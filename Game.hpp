#pragma once

#include "GL.hpp"

#include <SDL.h>
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

#include <vector>

// The 'Game' struct holds all of the game-relevant state,
// and is called by the main loop.

struct Game {
	//Game creates OpenGL resources (i.e. vertex buffer objects) in its
	//constructor and frees them in its destructor.
	Game();
	~Game();

	//handle_event is called when new mouse or keyboard events are received:
	// (note that this might be many times per frame or never)
	//The function should return 'true' if it handled the event.
	bool handle_event(SDL_Event const &evt, glm::uvec2 window_size);

	//update is called at the start of a new frame, after events are handled:
	void update(float elapsed);

	//draw is called after update:
	void draw(glm::uvec2 drawable_size);

	bool is_finished();

	void reset();

	void reset_newboard();

	//------- opengl resources -------

	//shader program that draws lit objects with vertex colors:
	struct {
		GLuint program = -1U; //program object

		//uniform locations:
		GLuint object_to_clip_mat4 = -1U;
		GLuint object_to_light_mat4x3 = -1U;
		GLuint normal_to_light_mat3 = -1U;
		GLuint sun_direction_vec3 = -1U;
		GLuint sun_color_vec3 = -1U;
		GLuint sky_direction_vec3 = -1U;
		GLuint sky_color_vec3 = -1U;

		//attribute locations:
		GLuint Position_vec4 = -1U;
		GLuint Normal_vec3 = -1U;
		GLuint Color_vec4 = -1U;
	} simple_shading;

	//mesh data, stored in a vertex buffer:
	GLuint meshes_vbo = -1U; //vertex buffer holding mesh data

	//The location of each mesh in the meshes vertex buffer:
	struct Mesh {
		GLint first = 0;
		GLsizei count = 0;
	};

	Mesh tile_mesh;
	Mesh cursor_mesh;
	Mesh doll_mesh;
	Mesh egg_mesh;
	Mesh cube_mesh;
	Mesh darkpiece_mesh;
	Mesh lightpiece_mesh;
	Mesh selectedtile_mesh;
	Mesh gameovertext_mesh;
	Mesh starttext_mesh;
	Mesh currentscoretext_mesh;
	Mesh bestscoretext_mesh;

	Mesh num0_mesh;
	Mesh num1_mesh;
	Mesh num2_mesh;
	Mesh num3_mesh;
	Mesh num4_mesh;
	Mesh num5_mesh;
	Mesh num6_mesh;
	Mesh num7_mesh;
	Mesh num8_mesh;
	Mesh num9_mesh;

	void merge_row_right(int row);
	void merge_row_left(int row);
	void merge_col_down(int col);
	void merge_col_up(int col);

	Mesh get_first_digit(uint32_t score);
	Mesh get_second_digit(uint32_t score);

	bool vertical_direction;
	bool game_over;
	bool show_start_screen;

	uint32_t lowscore;
	uint32_t currentscore;
	uint32_t seed;

	GLuint meshes_for_simple_shading_vao = -1U; //vertex array object that describes how to connect the meshes_vbo to the simple_shading_program

	//------- game state -------

	glm::uvec2 board_size = glm::uvec2(4,4);
	std::vector< Mesh const * > board_meshes;
	std::vector< glm::quat > board_rotations;

	glm::uvec2 cursor = glm::vec2(0,0);

	struct {
		bool roll_left = false;
		bool roll_right = false;
		bool roll_up = false;
		bool roll_down = false;
	} controls;

};
