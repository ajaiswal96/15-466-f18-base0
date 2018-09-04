#include "Game.hpp"

#include "gl_errors.hpp" //helper for dumpping OpenGL error messages
#include "read_chunk.hpp" //helper for reading a vector of structures from a file
#include "data_path.hpp" //helper to get paths relative to executable

#include <glm/gtc/type_ptr.hpp>

#include <iostream>
#include <fstream>
#include <map>
#include <cstddef>
#include <random>
#include <queue>

//helper defined later; throws if shader compilation fails:
static GLuint compile_shader(GLenum type, std::string const &source);

Game::Game() {
	{ //create an opengl program to perform sun/sky (well, directional+hemispherical) lighting:
		GLuint vertex_shader = compile_shader(GL_VERTEX_SHADER,
			"#version 330\n"
			"uniform mat4 object_to_clip;\n"
			"uniform mat4x3 object_to_light;\n"
			"uniform mat3 normal_to_light;\n"
			"layout(location=0) in vec4 Position;\n" //note: layout keyword used to make sure that the location-0 attribute is always bound to something
			"in vec3 Normal;\n"
			"in vec4 Color;\n"
			"out vec3 position;\n"
			"out vec3 normal;\n"
			"out vec4 color;\n"
			"void main() {\n"
			"	gl_Position = object_to_clip * Position;\n"
			"	position = object_to_light * Position;\n"
			"	normal = normal_to_light * Normal;\n"
			"	color = Color;\n"
			"}\n"
		);

		GLuint fragment_shader = compile_shader(GL_FRAGMENT_SHADER,
			"#version 330\n"
			"uniform vec3 sun_direction;\n"
			"uniform vec3 sun_color;\n"
			"uniform vec3 sky_direction;\n"
			"uniform vec3 sky_color;\n"
			"in vec3 position;\n"
			"in vec3 normal;\n"
			"in vec4 color;\n"
			"out vec4 fragColor;\n"
			"void main() {\n"
			"	vec3 total_light = vec3(0.0, 0.0, 0.0);\n"
			"	vec3 n = normalize(normal);\n"
			"	{ //sky (hemisphere) light:\n"
			"		vec3 l = sky_direction;\n"
			"		float nl = 0.5 + 0.5 * dot(n,l);\n"
			"		total_light += nl * sky_color;\n"
			"	}\n"
			"	{ //sun (directional) light:\n"
			"		vec3 l = sun_direction;\n"
			"		float nl = max(0.0, dot(n,l));\n"
			"		total_light += nl * sun_color;\n"
			"	}\n"
			"	fragColor = vec4(color.rgb * total_light, color.a);\n"
			"}\n"
		);

		//set the selected direction to be vertical by default
		vertical_direction = true;

		simple_shading.program = glCreateProgram();
		glAttachShader(simple_shading.program, vertex_shader);
		glAttachShader(simple_shading.program, fragment_shader);
		//shaders are reference counted so this makes sure they are freed after program is deleted:
		glDeleteShader(vertex_shader);
		glDeleteShader(fragment_shader);

		//link the shader program and throw errors if linking fails:
		glLinkProgram(simple_shading.program);
		GLint link_status = GL_FALSE;
		glGetProgramiv(simple_shading.program, GL_LINK_STATUS, &link_status);
		if (link_status != GL_TRUE) {
			std::cerr << "Failed to link shader program." << std::endl;
			GLint info_log_length = 0;
			glGetProgramiv(simple_shading.program, GL_INFO_LOG_LENGTH, &info_log_length);
			std::vector< GLchar > info_log(info_log_length, 0);
			GLsizei length = 0;
			glGetProgramInfoLog(simple_shading.program, GLsizei(info_log.size()), &length, &info_log[0]);
			std::cerr << "Info log: " << std::string(info_log.begin(), info_log.begin() + length);
			throw std::runtime_error("failed to link program");
		}
	}

	{ //read back uniform and attribute locations from the shader program:
		simple_shading.object_to_clip_mat4 = glGetUniformLocation(simple_shading.program, "object_to_clip");
		simple_shading.object_to_light_mat4x3 = glGetUniformLocation(simple_shading.program, "object_to_light");
		simple_shading.normal_to_light_mat3 = glGetUniformLocation(simple_shading.program, "normal_to_light");

		simple_shading.sun_direction_vec3 = glGetUniformLocation(simple_shading.program, "sun_direction");
		simple_shading.sun_color_vec3 = glGetUniformLocation(simple_shading.program, "sun_color");
		simple_shading.sky_direction_vec3 = glGetUniformLocation(simple_shading.program, "sky_direction");
		simple_shading.sky_color_vec3 = glGetUniformLocation(simple_shading.program, "sky_color");

		simple_shading.Position_vec4 = glGetAttribLocation(simple_shading.program, "Position");
		simple_shading.Normal_vec3 = glGetAttribLocation(simple_shading.program, "Normal");
		simple_shading.Color_vec4 = glGetAttribLocation(simple_shading.program, "Color");
	}

	struct Vertex {
		glm::vec3 Position;
		glm::vec3 Normal;
		glm::u8vec4 Color;
	};
	static_assert(sizeof(Vertex) == 28, "Vertex should be packed.");

	{ //load mesh data from a binary blob:
		std::ifstream blob(data_path("meshes.blob"), std::ios::binary);
		//The blob will be made up of three chunks:
		// the first chunk will be vertex data (interleaved position/normal/color)
		// the second chunk will be characters
		// the third chunk will be an index, mapping a name (range of characters) to a mesh (range of vertex data)

		//read vertex data:
		std::vector< Vertex > vertices;
		read_chunk(blob, "dat0", &vertices);

		//read character data (for names):
		std::vector< char > names;
		read_chunk(blob, "str0", &names);

		//read index:
		struct IndexEntry {
			uint32_t name_begin;
			uint32_t name_end;
			uint32_t vertex_begin;
			uint32_t vertex_end;
		};
		static_assert(sizeof(IndexEntry) == 16, "IndexEntry should be packed.");

		std::vector< IndexEntry > index_entries;
		read_chunk(blob, "idx0", &index_entries);

		if (blob.peek() != EOF) {
			std::cerr << "WARNING: trailing data in meshes file." << std::endl;
		}

		//upload vertex data to the graphics card:
		glGenBuffers(1, &meshes_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, meshes_vbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex) * vertices.size(), vertices.data(), GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		//create map to store index entries:
		std::map< std::string, Mesh > index;
		for (IndexEntry const &e : index_entries) {
			if (e.name_begin > e.name_end || e.name_end > names.size()) {
				throw std::runtime_error("invalid name indices in index.");
			}
			if (e.vertex_begin > e.vertex_end || e.vertex_end > vertices.size()) {
				throw std::runtime_error("invalid vertex indices in index.");
			}
			Mesh mesh;
			mesh.first = e.vertex_begin;
			mesh.count = e.vertex_end - e.vertex_begin;
			auto ret = index.insert(std::make_pair(
				std::string(names.begin() + e.name_begin, names.begin() + e.name_end),
				mesh));
			if (!ret.second) {
				throw std::runtime_error("duplicate name in index.");
			}
		}

		//look up into index map to extract meshes:
		auto lookup = [&index](std::string const &name) -> Mesh {
			auto f = index.find(name);
			if (f == index.end()) {
				throw std::runtime_error("Mesh named '" + name + "' does not appear in index.");
			}
			return f->second;
		};
		tile_mesh = lookup("Tile");
		cursor_mesh = lookup("Cursor");
		doll_mesh = lookup("Doll");
		egg_mesh = lookup("Egg");
		cube_mesh = lookup("Cube");
		darkpiece_mesh = lookup("DarkPiece");
		lightpiece_mesh = lookup("LightPiece");
		selectedtile_mesh = lookup("SelectedTile");
		gameovertext_mesh = lookup("GameOverText");
		starttext_mesh = lookup("StartText");
		currentscoretext_mesh = lookup("CurrentScoreText");
		bestscoretext_mesh = lookup("BestScoreText");

		num0_mesh = lookup("NUM0");
		num1_mesh = lookup("NUM1");
		num2_mesh = lookup("NUM2");
		num3_mesh = lookup("NUM3");
		num4_mesh = lookup("NUM4");
		num5_mesh = lookup("NUM5");
		num6_mesh = lookup("NUM6");
		num7_mesh = lookup("NUM7");
		num8_mesh = lookup("NUM8");
		num9_mesh = lookup("NUM9");
	}

	{ //create vertex array object to hold the map from the mesh vertex buffer to shader program attributes:
		glGenVertexArrays(1, &meshes_for_simple_shading_vao);
		glBindVertexArray(meshes_for_simple_shading_vao);
		glBindBuffer(GL_ARRAY_BUFFER, meshes_vbo);
		//note that I'm specifying a 3-vector for a 4-vector attribute here, and this is okay to do:
		glVertexAttribPointer(simple_shading.Position_vec4, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (GLbyte *)0 + offsetof(Vertex, Position));
		glEnableVertexAttribArray(simple_shading.Position_vec4);
		if (simple_shading.Normal_vec3 != -1U) {
			glVertexAttribPointer(simple_shading.Normal_vec3, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (GLbyte *)0 + offsetof(Vertex, Normal));
			glEnableVertexAttribArray(simple_shading.Normal_vec3);
		}
		if (simple_shading.Color_vec4 != -1U) {
			glVertexAttribPointer(simple_shading.Color_vec4, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(Vertex), (GLbyte *)0 + offsetof(Vertex, Color));
			glEnableVertexAttribArray(simple_shading.Color_vec4);
		}
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}

	GL_ERRORS();
	board_meshes.reserve(board_size.x * board_size.y);
	board_rotations.reserve(board_size.x * board_size.y);
	seed = 0xbead1234;
	lowscore = 0;
	show_start_screen = true;
	reset();
}

Game::~Game() {
	glDeleteVertexArrays(1, &meshes_for_simple_shading_vao);
	meshes_for_simple_shading_vao = -1U;

	glDeleteBuffers(1, &meshes_vbo);
	meshes_vbo = -1U;

	glDeleteProgram(simple_shading.program);
	simple_shading.program = -1U;

	GL_ERRORS();
}

bool Game::handle_event(SDL_Event const &evt, glm::uvec2 window_size) {
	//ignore any keys that are the result of automatic key repeat:
	if (evt.type == SDL_KEYDOWN && evt.key.repeat) {
		return false;
	}

	if (show_start_screen) {
		if (evt.type == SDL_KEYDOWN || evt.type == SDL_KEYUP) {
			if (evt.key.keysym.scancode == SDL_SCANCODE_SPACE) {
				show_start_screen = false;
				return true;
			}
		}
		return false;
	}
	//handle tracking the state of WSAD for roll control:
	if (evt.type == SDL_KEYDOWN || evt.type == SDL_KEYUP) {
		if (evt.key.keysym.scancode == SDL_SCANCODE_W) {
			controls.roll_up = (evt.type == SDL_KEYDOWN);
			if (evt.type == SDL_KEYDOWN && vertical_direction) {
				merge_col_up(cursor.x);
				currentscore++;
			}
			return true;
		} else if (evt.key.keysym.scancode == SDL_SCANCODE_S) {
			controls.roll_down = (evt.type == SDL_KEYDOWN);
			if(evt.type == SDL_KEYDOWN && vertical_direction) {
				merge_col_down(cursor.x);
				currentscore++;
			}
			return true;
		} else if (evt.key.keysym.scancode == SDL_SCANCODE_A) {
			controls.roll_left = (evt.type == SDL_KEYDOWN);
			if(evt.type == SDL_KEYDOWN && !vertical_direction) {
				merge_row_left(cursor.y);
				currentscore++;
			}
			return true;
		} else if (evt.key.keysym.scancode == SDL_SCANCODE_D) {
			controls.roll_right = (evt.type == SDL_KEYDOWN);
			if(evt.type == SDL_KEYDOWN && !vertical_direction) {
				merge_row_right(cursor.y);
				currentscore++;
			}
			return true;
		} else if (evt.key.keysym.scancode == SDL_SCANCODE_R) {
			if (evt.type == SDL_KEYDOWN){
				reset();
			}
			return true;
		}
	}
	//move cursor on L/R/U/D press:
	if (evt.type == SDL_KEYDOWN && evt.key.repeat == 0) {
		if (evt.key.keysym.scancode == SDL_SCANCODE_LEFT) {
			vertical_direction = true;
			if (cursor.x > 0) {
				cursor.x -= 1;
			}
			return true;
		} else if (evt.key.keysym.scancode == SDL_SCANCODE_RIGHT) {
			vertical_direction = true;
			if (cursor.x + 1 < board_size.x) {
				cursor.x += 1;
			}
			return true;
		} else if (evt.key.keysym.scancode == SDL_SCANCODE_UP) {
			vertical_direction = false;
			if (cursor.y + 1 < board_size.y) {
				cursor.y += 1;
			}
			return true;
		} else if (evt.key.keysym.scancode == SDL_SCANCODE_DOWN) {
			vertical_direction = false;
			if (cursor.y > 0) {
				cursor.y -= 1;
			}
			return true;
		}
	}

	if (evt.type == SDL_KEYDOWN && evt.key.repeat == 0) {
		if (evt.key.keysym.scancode == SDL_SCANCODE_SPACE) {
			reset_newboard();
			return true;
		}
	}
	return false;
}

void Game::update(float elapsed) {
	//if the roll keys are pressed, rotate everything on the same row or column as the cursor:
	game_over = is_finished();
	if (game_over) {
		if (currentscore!=0) {
			if (lowscore == 0 || currentscore < lowscore) {
				lowscore = currentscore;
			}
		}
	}

	glm::quat dr = glm::quat(1.0f, 0.0f, 0.0f, 0.0f);
	float amt = elapsed * 1.0f;
	//
	// if (!vertical_direction) {
	// 	dr = glm::angleAxis(amt, glm::vec3(1.0f, 0.0f, 0.0f)) * dr;
	// 	if (dr != glm::quat()) {
	// 		for (uint32_t x = 0; x < board_size.x; ++x) {
	// 			uint32_t index = cursor.x*board_size.y + x;
	// 			std::cout<<"index: "<<index<<"\n";
	// 			glm::quat &r = board_rotations[index];
	// 			r = glm::normalize(dr * r);
	// 		}
	// 	}
	// }

	if (controls.roll_left) {
		dr = glm::angleAxis(amt, glm::vec3(0.0f, 1.0f, 0.0f)) * dr;
	}
	if (controls.roll_right) {
		dr = glm::angleAxis(-amt, glm::vec3(0.0f, 1.0f, 0.0f)) * dr;
	}
	if (controls.roll_up) {
		dr = glm::angleAxis(amt, glm::vec3(1.0f, 0.0f, 0.0f)) * dr;
	}
	if (controls.roll_down) {
		dr = glm::angleAxis(-amt, glm::vec3(1.0f, 0.0f, 0.0f)) * dr;
	}
	if (dr != glm::quat()) {
		for (uint32_t x = 0; x < board_size.x; ++x) {
			glm::quat &r = board_rotations[cursor.y * board_size.x + x];
			r = glm::normalize(dr * r);
		}
		for (uint32_t y = 0; y < board_size.y; ++y) {
			if (y != cursor.y) {
				glm::quat &r = board_rotations[y * board_size.x + cursor.x];
				r = glm::normalize(dr * r);
			}
		}
	}
}

void Game::draw(glm::uvec2 drawable_size) {
	//Set up a transformation matrix to fit the board in the window:
	glm::mat4 world_to_clip;
	{
		float aspect = float(drawable_size.x) / float(drawable_size.y);

		//want scale such that board * scale fits in [-aspect,aspect]x[-1.0,1.0] screen box:
		float scale = glm::min(
			2.0f * aspect / float(board_size.x),
			2.0f / float(board_size.y)
		);

		//center of board will be placed at center of screen:
		glm::vec2 center = 0.5f * glm::vec2(board_size);

		//NOTE: glm matrices are specified in column-major order
		world_to_clip = glm::mat4(
			scale / aspect, 0.0f, 0.0f, 0.0f,
			0.0f, scale, 0.0f, 0.0f,
			0.0f, 0.0f,-1.0f, 0.0f,
			-(scale / aspect) * center.x, -scale * center.y, 0.0f, 1.0f
		);
	}

	//set up graphics pipeline to use data from the meshes and the simple shading program:
	glBindVertexArray(meshes_for_simple_shading_vao);
	glUseProgram(simple_shading.program);

	glUniform3fv(simple_shading.sun_color_vec3, 1, glm::value_ptr(glm::vec3(0.81f, 0.81f, 0.76f)));
	glUniform3fv(simple_shading.sun_direction_vec3, 1, glm::value_ptr(glm::normalize(glm::vec3(-0.2f, 0.2f, 1.0f))));
	glUniform3fv(simple_shading.sky_color_vec3, 1, glm::value_ptr(glm::vec3(0.2f, 0.2f, 0.3f)));
	glUniform3fv(simple_shading.sky_direction_vec3, 1, glm::value_ptr(glm::vec3(0.0f, 1.0f, 0.0f)));

	//helper function to draw a given mesh with a given transformation:
	auto draw_mesh = [&](Mesh const &mesh, glm::mat4 const &object_to_world) {
		//set up the matrix uniforms:
		if (simple_shading.object_to_clip_mat4 != -1U) {
			glm::mat4 object_to_clip = world_to_clip * object_to_world;
			glUniformMatrix4fv(simple_shading.object_to_clip_mat4, 1, GL_FALSE, glm::value_ptr(object_to_clip));
		}
		if (simple_shading.object_to_light_mat4x3 != -1U) {
			glUniformMatrix4x3fv(simple_shading.object_to_light_mat4x3, 1, GL_FALSE, glm::value_ptr(object_to_world));
		}
		if (simple_shading.normal_to_light_mat3 != -1U) {
			//NOTE: if there isn't any non-uniform scaling in the object_to_world matrix, then the inverse transpose is the matrix itself, and computing it wastes some CPU time:
			glm::mat3 normal_to_world = glm::inverse(glm::transpose(glm::mat3(object_to_world)));
			glUniformMatrix3fv(simple_shading.normal_to_light_mat3, 1, GL_FALSE, glm::value_ptr(normal_to_world));
		}

		//draw the mesh:
		glDrawArrays(GL_TRIANGLES, mesh.first, mesh.count);
	};

	for (uint32_t y = 0; y < board_size.y; ++y) {
		for (uint32_t x = 0; x < board_size.x; ++x) {
			draw_mesh(tile_mesh,
				glm::mat4(
					1.0f, 0.0f, 0.0f, 0.0f,
					0.0f, 1.0f, 0.0f, 0.0f,
					0.0f, 0.0f, 1.0f, 0.0f,
					x+0.5f, y+0.5f,-0.5f, 1.0f
				)
			);
			if (vertical_direction) {
				if(!game_over && !show_start_screen) {
					draw_mesh(selectedtile_mesh,
						glm::mat4(
							1.0f, 0.0f, 0.0f, 0.0f,
							0.0f, 1.0f, 0.0f, 0.0f,
							0.0f, 0.0f, 1.0f, 0.0f,
							cursor.x+0.5f, y+0.5f, 0.0f, 1.0f
						)
					);
				}
			}
			else {
				if (!game_over && !show_start_screen) {
					draw_mesh(selectedtile_mesh,
						glm::mat4(
							1.0f, 0.0f, 0.0f, 0.0f,
							0.0f, 1.0f, 0.0f, 0.0f,
							0.0f, 0.0f, 1.0f, 0.0f,
							x+0.5f, cursor.y+0.5f, 0.0f, 1.0f
						)
					);
				}

			}

			const Mesh* gamepiece_to_draw = board_meshes[y*board_size.x+x];
			if(gamepiece_to_draw) {
				draw_mesh(*gamepiece_to_draw,
						glm::mat4(
							1.0f, 0.0f, 0.0f, 0.0f,
							0.0f, 1.0f, 0.0f, 0.0f,
							0.0f, 0.0f, 1.0f, 0.0f,
							x+0.5f, y+0.5f, 0.0f, 1.0f
						)
						* glm::mat4_cast(board_rotations[y*board_size.x+x])
					);
				}
			}
		}

		if (show_start_screen) {
			draw_mesh(starttext_mesh,
					glm::mat4(
						1.0f, 0.0f, 0.0f, 0.0f,
						0.0f, 1.0f, 0.0f, 0.0f,
						0.0f, 0.0f, 1.0f, 0.0f,
						float(board_size.x-1)/2 - 0.5f, float(board_size.y-1)/2 + 0.5f, 1.0f, 1.0f
					)
				);
		}



		if (game_over) {
			draw_mesh(gameovertext_mesh,
				glm::mat4(
					1.0f, 0.0f, 0.0f, 0.0f,
					0.0f, 1.0f, 0.0f, 0.0f,
					0.0f, 0.0f, 1.0f, 0.0f,
					float(board_size.x-1)/2, float(board_size.y-1)/2, 0.1f, 1.0f
				)
			);
		}

		Mesh currscore_firstdigit_mesh = get_first_digit(currentscore);
		Mesh currscore_seconddigit_mesh = get_second_digit(currentscore);
		Mesh lowscore_firstdigit_mesh = get_first_digit(lowscore);
		Mesh lowscore_seconddigit_mesh = get_second_digit(lowscore);

		draw_mesh(currentscoretext_mesh,
			glm::mat4(
				1.0f, 0.0f, 0.0f, 0.0f,
				0.0f, 1.0f, 0.0f, 0.0f,
				0.0f, 0.0f, 1.0f, 0.0f,
				float(board_size.x) + 0.1f, float(board_size.y-0.7f), 0.0f, 1.0f
			)
		);

		draw_mesh(currscore_firstdigit_mesh,
			glm::mat4(
				1.0f, 0.0f, 0.0f, 0.0f,
				0.0f, 1.0f, 0.0f, 0.0f,
				0.0f, 0.0f, 1.0f, 0.0f,
				float(board_size.x) + 0.1f, float(board_size.y-1), 0.0f, 1.0f
			)
		);
		draw_mesh(currscore_seconddigit_mesh,
			glm::mat4(
				1.0f, 0.0f, 0.0f, 0.0f,
				0.0f, 1.0f, 0.0f, 0.0f,
				0.0f, 0.0f, 1.0f, 0.0f,
				float(board_size.x) + 0.3f, float(board_size.y-1), 0.0f, 1.0f
			)
		);

		draw_mesh(bestscoretext_mesh,
			glm::mat4(
				1.0f, 0.0f, 0.0f, 0.0f,
				0.0f, 1.0f, 0.0f, 0.0f,
				0.0f, 0.0f, 1.0f, 0.0f,
				float(board_size.x) + 0.1f, float(board_size.y-1.7), 0.0f, 1.0f
			)
		);

		draw_mesh(lowscore_firstdigit_mesh,
			glm::mat4(
				1.0f, 0.0f, 0.0f, 0.0f,
				0.0f, 1.0f, 0.0f, 0.0f,
				0.0f, 0.0f, 1.0f, 0.0f,
				float(board_size.x) + 0.1f, float(board_size.y-2), 0.0f, 1.0f
			)
		);
		draw_mesh(lowscore_seconddigit_mesh,
			glm::mat4(
				1.0f, 0.0f, 0.0f, 0.0f,
				0.0f, 1.0f, 0.0f, 0.0f,
				0.0f, 0.0f, 1.0f, 0.0f,
				float(board_size.x) + 0.3f, float(board_size.y-2), 0.0f, 1.0f
			)
		);

	glUseProgram(0);

	GL_ERRORS();
}


void Game::merge_row_left(int row){
	bool element_deleted = false;
	uint32_t first_index = row*board_size.y;
	const Mesh* first_mesh = board_meshes[first_index];
	for (uint32_t x = 1; x < board_size.x; ++x) {
		uint32_t curr_index = row*board_size.y + x;
		const Mesh* curr_mesh = board_meshes[curr_index];
		if (first_mesh == curr_mesh) {
				board_meshes[curr_index] = nullptr;
				board_meshes[first_index] = nullptr;
				element_deleted = true;
			} else {
				break;
			}
	}
	if (element_deleted) {
		uint32_t num_deleted = 0;
		for (uint32_t x = 0; x < board_size.x; ++x) {
			uint32_t current_index = row*board_size.y + x;
			if (!board_meshes[current_index]){
				num_deleted ++;
			}
			else {
				board_meshes[current_index - num_deleted] = board_meshes[current_index];
				board_meshes[current_index] = nullptr;
			}

		}
	}

	std::queue<const Mesh*> non_null_element_queue;
	for (uint32_t x =0; x < board_size.x; ++x) {
		uint32_t current_index = row*board_size.y + x;
		if (board_meshes[current_index]) {
			non_null_element_queue.push(board_meshes[current_index]);
		}
	}
	for (uint32_t x =0; x < board_size.x; ++x) {
		uint32_t current_index = row*board_size.y + x;
		if (!non_null_element_queue.empty()) {
			board_meshes[current_index] = non_null_element_queue.front();
			non_null_element_queue.pop();
		} else {
			board_meshes[current_index] = nullptr;
		}
	}
}

void Game::merge_row_right(int row){
	bool element_deleted = false;
	uint32_t first_index = row*board_size.y + (board_size.x-1);
	const Mesh* first_mesh = board_meshes[first_index];
	for (int x = board_size.x-2; x >=0; --x) {
		uint32_t curr_index = row*board_size.y + x;
		const Mesh* curr_mesh = board_meshes[curr_index];
		if (first_mesh == curr_mesh) {
			board_meshes[curr_index] = nullptr;
			board_meshes[first_index] = nullptr;
			element_deleted = true;
		} else {
			break;
		}
	}
	if (element_deleted){
		uint32_t num_deleted = 0;
		for (int x = board_size.x-1; x>=0; --x) {
			uint32_t current_index = row*board_size.y + x;
			if (!board_meshes[current_index]){
				num_deleted ++;
			}
			else {
				board_meshes[current_index + num_deleted] = board_meshes[current_index];
				board_meshes[current_index] = nullptr;
			}
		}
	}

	std::queue<const Mesh*> non_null_element_queue;
	for (int x = board_size.x-1; x>=0; --x) {
		uint32_t current_index = row*board_size.y + x;
		if (board_meshes[current_index]) {
			non_null_element_queue.push(board_meshes[current_index]);
		}
	}
	for (int x = board_size.x-1; x>=0; --x) {
		uint32_t current_index = row*board_size.y + x;
		if (!non_null_element_queue.empty()){
			board_meshes[current_index] = non_null_element_queue.front();
			non_null_element_queue.pop();
		} else {
			board_meshes[current_index] = nullptr;
		}
	}
}

void Game::merge_col_down(int col){
	bool element_deleted = false;
	uint32_t first_index = col;
	const Mesh* first_mesh = board_meshes[first_index];
	for (uint32_t y = 1; y < board_size.y; ++y) {
		uint32_t curr_index = y*board_size.x+col;
		const Mesh* curr_mesh = board_meshes[curr_index];
		if (first_mesh == curr_mesh) {
			board_meshes[first_index] = nullptr;
			board_meshes[curr_index] = nullptr;
			element_deleted = true;
		} else {
			break;
		}
	}
	if (element_deleted){
		uint32_t num_deleted = 0;
		for (uint32_t y = 0; y < board_size.y; ++y){
			uint32_t current_index = y*board_size.x+col;
			if (!board_meshes[current_index]){
				num_deleted ++;
			}
			else {
				board_meshes[current_index - board_size.x*num_deleted] = board_meshes[current_index];
				board_meshes[current_index] = nullptr;
			}
		}
	}

	std::queue<const Mesh*> non_null_element_queue;
	for (uint32_t y = 0; y < board_size.y; ++y){
		uint32_t current_index = y*board_size.x+col;
		if (board_meshes[current_index]) {
			non_null_element_queue.push(board_meshes[current_index]);
		}
	}
	for (uint32_t y = 0; y < board_size.y; ++y){
		uint32_t current_index = y*board_size.x+col;
		if(!non_null_element_queue.empty()){
			board_meshes[current_index] = non_null_element_queue.front();
			non_null_element_queue.pop();
		} else {
			board_meshes[current_index] = nullptr;
		}
	}
}

void Game::merge_col_up(int col){
	bool element_deleted = false;
	uint32_t first_index = (board_size.y -1)*board_size.x + col;
	const Mesh* first_mesh = board_meshes[first_index];
	for (int y = board_size.y-2; y >=0; --y) {
		uint32_t curr_index = y*board_size.x+col;
		const Mesh* curr_mesh = board_meshes[curr_index];
		if (first_mesh == curr_mesh) {
			board_meshes[first_index] = nullptr;
			board_meshes[curr_index] = nullptr;
			element_deleted = true;
		} else {
			break;
		}
	}
	if (element_deleted){
		uint32_t num_deleted = 0;
		for (int y = board_size.y-1; y>=0; --y){
			uint32_t current_index = y*board_size.x+col;
			if (!board_meshes[current_index]){
				num_deleted ++;
			}
			else {
				board_meshes[current_index + board_size.x*num_deleted] = board_meshes[current_index];
				board_meshes[current_index] = nullptr;
			}
		}
	}

	std::queue<const Mesh*> non_null_element_queue;
	for (int y = board_size.y-1; y>=0; --y){
		uint32_t current_index = y*board_size.x+col;
		if (board_meshes[current_index]) {
			non_null_element_queue.push(board_meshes[current_index]);
		}
	}
	for (int y = board_size.y-1; y>=0; --y){
		uint32_t current_index = y*board_size.x+col;
		if(!non_null_element_queue.empty()){
			board_meshes[current_index] = non_null_element_queue.front();
			non_null_element_queue.pop();
		} else {
			board_meshes[current_index] = nullptr;
		}
	}
}

void Game::reset_newboard(){
	lowscore = 0;
	currentscore = 0;
	game_over = false;

	board_meshes.clear();
	board_rotations.clear();
	//get seed based on time for a new board
	//THIS LINE IS TAKEN FROM THE CPLUSPLUS REFERENCE WEBSITE: http://www.cplusplus.com/reference/random/mersenne_twister_engine/mersenne_twister_engine/
	seed = std::chrono::system_clock::now().time_since_epoch().count();
	reset();
}

void Game::reset(){
	currentscore = 0;
	game_over = false;

	board_meshes.clear();
	board_rotations.clear();
	//----------------
	//set up game board with meshes and rolls:

	std::mt19937 mt(seed);

	// std::vector< Mesh const * > meshes{ &doll_mesh, &egg_mesh, &cube_mesh, &darkpiece_mesh, &lightpiece_mesh};
	std::vector< Mesh const * > meshes{&darkpiece_mesh, &lightpiece_mesh};

	for (uint32_t i = 0; i < board_size.x * board_size.y; ++i) {
		board_meshes.emplace_back(meshes[mt()%meshes.size()]);
		board_rotations.emplace_back(glm::quat());
	}
}

bool Game::is_finished(){
	std::vector<Mesh const *> remaining_meshes;
	for (uint32_t i = 0; i < board_size.x * board_size.y; ++i) {
		if (board_meshes[i]) {
			remaining_meshes.push_back(board_meshes[i]);
		}
	}
	if (remaining_meshes.size() < 2){
		return true;
	} else if (remaining_meshes.size() == 2) {
		if(remaining_meshes[0] != remaining_meshes[1]) {
			return true;
		}
	}
	return false;
}

Game::Mesh Game::get_first_digit(uint32_t score) {
	uint32_t first_digit = score/10;
	if (first_digit > 9) {
		first_digit = 9;
	}
	switch(first_digit)
	{
		case 1: return num1_mesh;
		case 2: return num2_mesh;
		case 3: return num3_mesh;
		case 4: return num4_mesh;
		case 5: return num5_mesh;
		case 6: return num6_mesh;
		case 7: return num7_mesh;
		case 8: return num8_mesh;
		case 9: return num9_mesh;
		default: return num0_mesh;
	}
}

Game::Mesh Game::get_second_digit(uint32_t score) {
	uint32_t second_digit = score%10;
	switch(second_digit)
	{
		case 1: return num1_mesh;
		case 2: return num2_mesh;
		case 3: return num3_mesh;
		case 4: return num4_mesh;
		case 5: return num5_mesh;
		case 6: return num6_mesh;
		case 7: return num7_mesh;
		case 8: return num8_mesh;
		case 9: return num9_mesh;
		default: return num0_mesh;
	}
}


//create and return an OpenGL vertex shader from source:
static GLuint compile_shader(GLenum type, std::string const &source) {
	GLuint shader = glCreateShader(type);
	GLchar const *str = source.c_str();
	GLint length = GLint(source.size());
	glShaderSource(shader, 1, &str, &length);
	glCompileShader(shader);
	GLint compile_status = GL_FALSE;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &compile_status);
	if (compile_status != GL_TRUE) {
		std::cerr << "Failed to compile shader." << std::endl;
		GLint info_log_length = 0;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &info_log_length);
		std::vector< GLchar > info_log(info_log_length, 0);
		GLsizei length = 0;
		glGetShaderInfoLog(shader, GLsizei(info_log.size()), &length, &info_log[0]);
		std::cerr << "Info log: " << std::string(info_log.begin(), info_log.begin() + length);
		glDeleteShader(shader);
		throw std::runtime_error("Failed to compile shader.");
	}
	return shader;
}
