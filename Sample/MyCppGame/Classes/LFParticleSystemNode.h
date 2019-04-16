#ifndef __LiquidFun_EyeCandy__LFParticleSystemNode__
#define __LiquidFun_EyeCandy__LFParticleSystemNode__

#include "Box2D/Box2D.h"
#include "cocos2d.h"

class LFParticleSystemNode : public cocos2d::Node, public cocos2d::BlendProtocol
{
public:
    static LFParticleSystemNode* create(b2ParticleSystem* particleSystem, float ratio, const std::string&);

    virtual void draw(cocos2d::Renderer *renderer, const cocos2d::Mat4 &transform, uint32_t transformFlags) override;
    void onDraw(const cocos2d::Mat4 &transform, uint32_t transformFlags);

    const cocos2d::BlendFunc& getBlendFunc() const override;
    void setBlendFunc(const cocos2d::BlendFunc &var) override;

protected:
    LFParticleSystemNode();
    ~LFParticleSystemNode();

    bool init(b2ParticleSystem* particleSystem, float ratio, const std::string&);
    void setup();

    cocos2d::CustomCommand _customCommand;
    b2ParticleSystem* _particleSystem; // weak ref
    cocos2d::Mat4 _ratioTransform;
    float _ratio;

    GLfloat *_sizes;

    cocos2d::BlendFunc _blendFunc;
    cocos2d::Texture2D *_texture;
    cocos2d::RenderTexture *_renderTexture;
};

#endif /* defined(__LiquidFun_EyeCandy__LFParticleSystemNode__) */