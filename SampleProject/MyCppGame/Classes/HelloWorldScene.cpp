/****************************************************************************
 Copyright (c) 2017-2018 Xiamen Yaji Software Co., Ltd.
 
 http://www.cocos2d-x.org
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 ****************************************************************************/

#include "HelloWorldScene.h"
#include "SimpleAudioEngine.h"
#include "LFParticleSystemNode.h"
#include "Box2D/Box2D.h"
#include <cmath>
#include <ctime>

USING_NS_CC;

Scene* HelloWorld::createScene()
{
    return HelloWorld::create();
}

// Print useful error message instead of segfaulting when files are not there.
static void problemLoading(const char* filename)
{
    printf("Error while loading: %s\n", filename);
    printf("Depending on how you compiled you might have to add 'Resources/' in front of filenames in HelloWorldScene.cpp\n");
}

// on "init" you need to initialize your instance
bool HelloWorld::init()
{
    //////////////////////////////
    // 1. super init first
    if ( !Scene::init() )
    {
        return false;
    }

	world.reset(new b2World(b2Vec2(0, 0)));

	b2ParticleSystemDef def;
	def.destroyByAge = true;
	def.ejectionStrength = .1;
	def.powderStrength = .8;
	def.radius = 5;

	auto ps = world->CreateParticleSystem(&def);
	ps->SetMaxParticleCount(100);
	srand(time(nullptr));
    auto visibleSize = Director::getInstance()->getVisibleSize();

	schedule([ps, visibleSize](float dt) {
	   	// Create one particle per delta time

		b2ParticleDef def;

		auto p_x = rand() % (int)visibleSize.width;
		auto p_y = rand() % (int)visibleSize.height;
		def.position = b2Vec2(p_x, p_y);
		def.flags = b2ParticleFlag::b2_powderParticle;
		auto v_x = rand() % 5;
		auto v_y = rand() % 5;
		def.velocity = b2Vec2(v_x, v_y);
		def.lifetime = .5;

		b2ParticleColor c;
		c.a = 255;
		c.r = 255;
        c.g = 255;
        c.b = 255;

		def.color = c;

		ps->CreateParticle(def);

	}, "particleCreator");

	auto lf = LFParticleSystemNode::create(ps, 1.0f, "");
	addChild(lf);
    lf->setPosition(0, 0);

	scheduleUpdate();
    return true;
}

void HelloWorld::update(float dt) {
    // using the default particle iterations
	world->Step(dt, 8, 3);
}
