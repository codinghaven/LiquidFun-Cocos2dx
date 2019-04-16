Welcome to LiquidFun and Cocos2dx.

If you are looking for easy integration for LiquidFun and Cocos2dx than you have come to the right place.

This project was made to simplify the process of adding LiquidFun/Box2D into your project with cocos2d.
Instead of building and linking LiquidFun into architectural dependant libraries, use a simple `.cpp` file
that contains the whole source code of LiquidFun.

This version of cocos2d is 3.17.1.

External 3rd party libs are already installed with Box2d modified, namely prebuilt libraries and all .cpp
file from LiquidFun is removed except for the liquid.cpp file which contains everything.

When creating a new project using `cocos new...` using this cocos2dx version you will get LiquidFun
already integrated.

If you don't want to use this forked version of cocos2d, you will have to manually disable and remove
all Box2D prebuilt folders and static linking in the related build files, failing to do so will result in numerous
build errors due to multiple definitions.

The `liquid.cpp` file can be found under the `cocos2dx/external/Box2D` 

Creating a project is the same as you would regularly do using `cocos new`

After that you will have an empty project that's ready to go :)

Note: 
Windows is not yet supported to work with `cocos new`.

Update:
Android and IOS fully supported.
