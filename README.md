Welcome to LiquidFun and Cocos2dx.

If you are looking for easy integration for LiquidFun and Cocos2dx than you have come to the right place.

This project was made to simplify the process of adding LiquidFun/Box2D into your project with cocos2d.
Instead of building and linking LiquidFun into architectural dependant libraries, use a simple `.hpp` file
that contains the whole source code of LiquidFun and include it wherever you need it.

This version of cocos2d is 3.17.1.

External 3rd party libs are already installed except for Box2d (obviously).

When creating a new project using `cocos new...` nothing whatsover related to Box2d will be creating
so there's no static linking to any Box2d liibrary as there is none in this forked version of cocos2dx.
Instead when creating a project you will find a `liquid.hpp` file under the `Classes` folder.

If you don't want to use this forked version of cocos2d, you will have to manually disable and remove
and Box2D configurations, you don't have to but it is advised.

Instructions:
First read instructions on cocos2d-x main page 
Clone the project and and run the regular `cocos new...` command.

After that you will have an empty project that's ready to go :)
