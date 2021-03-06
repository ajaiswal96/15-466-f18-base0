# Game Information
Title: Slidey Ball

Author: Anubhav Jaiswal

Design Document: This game was inspired by the game in the following document. However certain improvements were made to the specified game.
1) Instead of just merging adjacent balls in a powerful slide, balls merge like they do in 2048 --> adjacent ones matching the edge merge together.
2) Instead of having a separate slide and powerful slide, there is a single merge. A slide is automatically performed if balls cannot be combined. Otherwise a 'powerful slide' is performed, followed by a slide. I wanted to match the dynamics of 2048, and this makes the gameplay more enjoyable.
3) Instead of merging arbitrary adjacent balls, only those that are adjacent to and match the edge to which they are being merged get merged.
3) The player has a choice of replaying a level ('r') or trying a new level ('space'). If the player replays a level, the best score is retained. If the player starts a new random level, the best score is reset.
4) A player can highlight a row or column using the arrow keys. If a column is highlighted, then only up and down merges can be performed. If a row is selected, only left and right merges can be performed. These merges are controlled via 'WASD'

Here is the original design doc which inspired this game:[Sliding Ball](http://graphics.cs.cmu.edu/courses/15-466-f18/game0-designs/guanghay/)

Screen Shot:

![Screen Shot](slideyball_screenshot.png)

Difficulties Encountered: The main difficulties with this game was learning OpenGL, Blender, and touching up on C++. It has been years since I used both C++ and OpenGL, so I spent a considerable amount of time going through the code and reading OpenGL documentation / watching videos about it's pipeline and what happens at each stage. Learning Blender wasn't too bad though it was annoying to create meshes for all the text needed for scores and menus. The other difficulty I encountered (mentioned in the good code section) was trying to make moves in the game more intuitive. In hindsight, I should have just stuck to the design doc for the sake of time.

Good Code: I felt the code that involved merging was done fairly well. The way the code is written may not be the best, though the logic is pretty good. I had first written slides/powerful slides as specified by the original game document, though they weren't as intuitive to use. Thus I used an iterative procedure by writing different merge functions and had friends play the game on each in order to come up with what we thought was the most intuitive. The code in the 'merge_row/col_left/right/up/down' methods illustrates this iterative process, and while not as concise as possible, it is super easy to debug.

# Controls
-Use the arrow keys to select a row or column in which you want to merge balls (the selected row/col will be highlighted in green).
-Once a row/column is selected use 'a'/'d' to merge balls left or right in a row and 'w'/'s' to merge balls up or down in a column
-Merging occurs similar to 2048. For example if a right merge occurs, all elements matching and adjacent to the element at the right edge of that row will be removed. After that all elements will be shifted to the right to fill up the space. If nothing can be merged, then all elements shift right to fill any gaps in that row. This process is the same for all directions of merges.
-The game ends when no balls can be merged with each other.
-The goal of the game is to make the fewest merges as possible (have the lowest score).
-Your score increases by 1 for each merge you ATTEMPT.
-The best score remains for a level, and you can retry a level by pressing 'r'
-You can select a new level to play by pressing 'space', though this will reset your best score

# Using This Base Code

Before you dive into the code, it helps to understand the overall structure of this repository.
- Files you should read and/or edit:
    - ```main.cpp``` creates the game window and contains the main loop. You should read through this file to understand what it's doing, but you shouldn't need to change things (other than window title and size).
    - ```Game.*pp``` declaration+definition for the Game struct. These files will contain the bulk of your code changes.
    - ```meshes/export-meshes.py``` exports meshes from a .blend file into a format usable by our game runtime. You will need to edit this file to add vertex color export code.
    - ```Jamfile``` responsible for telling FTJam how to build the project. If you add any additional .cpp files or want to change the name of your runtime executable you will need to modify this.
    - ```.gitignore``` ignores the ```objs/``` directory and the generated executable file. You will need to change it if your executable name changes. (If you find yourself changing it to ignore, e.g., your editor's swap files you should probably, instead be investigating making this change in the global git configuration.)
- Files you probably should at least glance at because they are useful:
    - ```read_chunk.hpp``` contains a function that reads a vector of structures prefixed by a magic number. It's surprising how many simple file formats you can create that only require such a function to access.
    - ```data_path.*pp``` contains a helper function that allows you to specify paths relative to the executable (instead of the current working directory). Very useful when loading assets.
	- ```gl_errors.hpp``` contains a function that checks for opengl error conditions. Also, the helpful macro ```GL_ERRORS()``` which calls ```gl_errors()``` with the current file and line number.
- Files you probably don't need to read or edit:
    - ```GL.hpp``` includes OpenGL prototypes without the namespace pollution of (e.g.) SDL's OpenGL header. It makes use of ```glcorearb.h``` and ```gl_shims.*pp``` to make this happen.
    - ```make-gl-shims.py``` does what it says on the tin. Included in case you are curious. You won't need to run it.

## Asset Build Instructions

In order to generate the ```dist/meshes.blob``` file, tell blender to execute the ```meshes/export-meshes.py``` script:

```
blender --background --python meshes/export-meshes.py -- meshes/meshes.blend dist/meshes.blob
```

There is a Makefile in the ```meshes``` directory that will do this for you.

## Runtime Build Instructions

The runtime code has been set up to be built with [FT Jam](https://www.freetype.org/jam/).

### Getting Jam

For more information on Jam, see the [Jam Documentation](https://www.perforce.com/documentation/jam-documentation) page at Perforce, which includes both reference documentation and a getting started guide.

On unixish OSs, Jam is available from your package manager:
```
	brew install ftjam #on OSX
	apt get ftjam #on Debian-ish Linux
```

On Windows, you can get a binary [from sourceforge](https://sourceforge.net/projects/freetype/files/ftjam/2.5.2/ftjam-2.5.2-win32.zip/download),
and put it somewhere in your `%PATH%`.
(Possibly: also set the `JAM_TOOLSET` variable to `VISUALC`.)

### Libraries

This code uses the [libSDL](https://www.libsdl.org/) library to create an OpenGL context, and the [glm](https://glm.g-truc.net) library for OpenGL-friendly matrix/vector types.
On MacOS and Linux, the code should work out-of-the-box if if you have these installed through your package manager.

If you are compiling on Windows or don't want to install these libraries globally there are pre-built library packages available in the
[kit-libs-linux](https://github.com/ixchow/kit-libs-linux),
[kit-libs-osx](https://github.com/ixchow/kit-libs-osx),
and [kit-libs-win](https://github.com/ixchow/kit-libs-win) repositories.
Simply clone into a subfolder and the build should work.

### Building

Open a terminal (or ```x64 Native Tools Command Prompt for VS 2017``` on Windows), change to the directory containing this code, and type:

```
jam
```

That's it. You can use ```jam -jN``` to run ```N``` parallel jobs if you'd like; ```jam -q``` to instruct jam to quit after the first error; ```jam -dx``` to show commands being executed; or ```jam main.o``` to build a specific file (in this case, main.cpp).  ```jam -h``` will print help on additional options.
