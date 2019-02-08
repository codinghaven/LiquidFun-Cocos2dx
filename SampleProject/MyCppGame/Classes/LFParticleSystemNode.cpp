#include "LFParticleSystemNode.h"

//
// --------- Custom Shaders ---------
//          Particle System
//
#define STRINGIFY(A)  #A

static const GLchar* _particleShaderVert = STRINGIFY(
attribute vec4 a_position;
uniform float u_size;

void main()
{
    gl_Position = CC_PMatrix * CC_MVMatrix * a_position;
    gl_PointSize = CC_MVMatrix[0][0] * u_size * 1.5;
}
);

// Fragment shader
static const GLchar* _particleShaderFrag = STRINGIFY(

#ifdef GL_ES
precision lowp float;
#endif

void main()
{
    gl_FragColor = texture2D(CC_Texture0, gl_PointCoord);
}
);

//
// --------- Custom Shaders ---------
//           Render Texture
//

// Vertex shader
static const GLchar* _renderTextureShaderVert = STRINGIFY(
attribute vec4 a_position;
attribute vec2 a_texCoord;
attribute vec4 a_color;

#ifdef GL_ES
varying lowp vec4 v_fragmentColor;
varying mediump vec2 v_texCoord;
#else
varying vec4 v_fragmentColor;
varying vec2 v_texCoord;
#endif

void main()
{
    gl_Position = CC_PMatrix * a_position;
    v_fragmentColor = a_color;
    v_texCoord = a_texCoord;
}
);

// Fragment shader
static const GLchar* _renderTextureShaderFrag = STRINGIFY(

#ifdef GL_ES
precision lowp float;
#endif

varying vec4 v_fragmentColor;
varying vec2 v_texCoord;
uniform float u_threshold_discard;
uniform float u_threshold_border;

void main()
{
    //    gl_FragColor = v_fragmentColor * texture2D(CC_Texture0, v_texCoord);
    vec4 color = v_fragmentColor * texture2D(CC_Texture0, v_texCoord);
    if( color.r < u_threshold_discard)
        // black or discard
        color = vec4(0,0,0,0);
    else if( color.r < u_threshold_border)
        // blue for the border
        color = vec4(0.2,0.2,0.9,1);
    else
        // white for the center
        color = vec4(1,1,1,1);
    gl_FragColor = color;
});



using namespace cocos2d;

LFParticleSystemNode::LFParticleSystemNode()
: _ratioTransform(Mat4::IDENTITY)
{
}

LFParticleSystemNode::~LFParticleSystemNode()
{
    free(_sizes);
}

LFParticleSystemNode* LFParticleSystemNode::create(b2ParticleSystem* particleSystem, float ratio, const std::string& texture_file)
{
    auto node = new LFParticleSystemNode;
    node->init(particleSystem, ratio, texture_file);
    node->autorelease();

    return node;
}

bool LFParticleSystemNode::init(b2ParticleSystem* particleSystem, float ratio, const std::string& texture_file)
{
    _particleSystem = particleSystem;
    _ratioTransform.scale(ratio, ratio, 1);
    _ratio = ratio;

    // own shader

    auto glprogram = GLProgram::createWithByteArrays(_particleShaderVert, _particleShaderFrag);

    GLProgramState *state = GLProgramState::getOrCreateWithGLProgram(glprogram);
    setGLProgramState(state);

    setup();


    _blendFunc = BlendFunc::ALPHA_NON_PREMULTIPLIED;
//    _blendFunc = BlendFunc::ADDITIVE;


    // create a render texture, this is what we are going to draw into
    auto s = Director::getInstance()->getWinSize();
    _renderTexture = cocos2d::RenderTexture::create(s.width, s.height, Texture2D::PixelFormat::RGBA8888);
    this->addChild(_renderTexture);
    _renderTexture->setAnchorPoint(Point::ANCHOR_MIDDLE);
    _renderTexture->setPosition(Point(s.width/2, s.height/2));

    // Change RenderTexture shader
    auto program = GLProgram::createWithByteArrays(_renderTextureShaderVert, _renderTextureShaderFrag);
    auto programState = GLProgramState::getOrCreateWithGLProgram(program);
    programState->setUniformFloat("u_threshold_discard", 0.15);
    programState->setUniformFloat("u_threshold_border", 0.3);

    _renderTexture->getSprite()->setGLProgramState(programState);

    return true;
}

void LFParticleSystemNode::setup()
{
    auto glProgramState = getGLProgramState();
    float size = _particleSystem->GetRadius()*2;

    // Update Positions
    glProgramState->setVertexAttribPointer("a_position", 2, GL_FLOAT, GL_FALSE, 0, _particleSystem->GetPositionBuffer());
    glProgramState->setUniformFloat("u_size", size);

    CHECK_GL_ERROR_DEBUG();
}

// draw

void LFParticleSystemNode::draw(Renderer *renderer, const Mat4 &transform, uint32_t transformFlags)
{
    _renderTexture->beginWithClear(0,0,0,0);

    _customCommand.init(_globalZOrder);
    _customCommand.func = CC_CALLBACK_0(LFParticleSystemNode::onDraw, this, transform, transformFlags);
    renderer->addCommand(&_customCommand);

    _renderTexture->end();
}

void LFParticleSystemNode::onDraw(const Mat4 &transform, uint32_t transformFlags)
{
    // transform everything to PTM_RATIO
    Mat4 newMV;
    Mat4::multiply(_modelViewTransform, _ratioTransform, &newMV);

    _glProgramState->apply(newMV);

    GL::blendFunc(_blendFunc.src, _blendFunc.dst);

    int totalParticles = _particleSystem->GetParticleCount();

    // Update Positions
    glDrawArrays(GL_POINTS, 0, totalParticles);
    CC_INCREMENT_GL_DRAWN_BATCHES_AND_VERTICES(1,totalParticles);

    CHECK_GL_ERROR_DEBUG();
}

// Blend Interface

/// blendFunc getter
const BlendFunc& LFParticleSystemNode::getBlendFunc() const
{
    return _blendFunc;
}
/// blendFunc setter
void LFParticleSystemNode::setBlendFunc(const BlendFunc &var)
{
    _blendFunc = var;
}
