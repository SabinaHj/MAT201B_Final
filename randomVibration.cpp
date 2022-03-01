// Karl Yerkes
// MAT201B
// 2022-01-04
// minimal app, ready for adapting..
//

#include "al/app/al_App.hpp"
#include "al/app/al_GUIDomain.hpp"
#include "al/graphics/al_Shapes.hpp"
#include "Gamma/Noise.h"
#include "Gamma/Delay.h"
#include "Gamma/Oscillator.h"
#include "Gamma/Filter.h"
#include "Gamma/Envelope.h"
#include "Gamma/Effects.h"
#include "Gamma/Analysis.h"
#include "al/math/al_Random.hpp"  // al::rnd::uniform(1, 2) , rnd::normal(..)

using namespace al;
using namespace gam;

struct MyApp : App {
  Parameter scale_x{"x", 1, 0, 2};
  Parameter scale_y{"y", 1, 0, 2};
  Parameter scale_z{"z", 1, 0, 2};
  Parameter gain{"gain", 0.5, 0, 1};

  Pan<> mPan;
  Sine<> mOsc1;
  Sine<> mOsc2;
  Sine<> mOsc3;
  Env<3> mAmpEnv;
  // envelope follower to connect audio output to graphics
  EnvFollow<> mEnvFollow;

  Mesh ball;
  NoisePink<> noise;
	Comb<float, ipl::Switchable> comb;
  double a = 0;
  double b = 0;
 	LFO<> mod;
	OnePole<> onePole;
  Vec3f rand_vec;
  
  void onCreate() override {
    addSphere(ball, 1, 100, 100);
    ball.decompress();
    ball.generateNormals();
 
    ball.color(al::rnd::uniform(),al::rnd::uniform(), al::rnd::uniform());
    
    nav().pos(0, 0, 5);
    comb.maxDelay(1./100);
  	comb.ipolType(ipl::ROUND);
  }
  void onAnimate(double dt) override {
    a += 0.29;
    b += 0.23;
    rand_vec = {scale_x.get() * al::rnd::uniform(0.5, 1.), scale_y.get() * al::rnd::uniform(0.5, 1.), scale_z.get() * al::rnd::uniform(0.5, 1.)};
  }

  void onDraw(Graphics& g) override {
    g.clear(0.2);
    g.depthTesting(true);
    g.lighting(true);
    g.rotate(a, Vec3f(0, 1, 0));
    g.rotate(b, Vec3f(1));
    g.scale(rand_vec);
//    g.scale( scale_x.get() * rand(), scale_y.get() * rand(),scale_z.get() * rand());
   g.meshColor();
  // g.color(0, 1, 1);
    g.draw(ball);
  }

  void onSound(AudioIOData &io) override {
    mOsc1.freq(rand_vec[0]*100+50);
    mOsc2.freq(rand_vec[1]*500+10);
    mOsc3.freq(rand_vec[2]*8000);

    // mOsc1.freq(scale_x.get()*100+50);
    // mOsc2.freq(scale_y.get()*1000+50);
    // mOsc3.freq(scale_z.get()*3000+50);


//		mod.period(scale_x.get());
//		comb.delay(mod.triU() * scale_y.get() /1000 + scale_z.get()/10);
    while (io()) {
      io.out(0) = mOsc1() + mOsc3() * gain;
      io.out(1) = mOsc1() + mOsc2() * (1-gain);

    }
  }
  void onInit() override {
    auto GUIdomain = GUIDomain::enableGUI(defaultWindowDomain());
    auto& gui = GUIdomain->newGUI();
    gui.add(scale_x);
    gui.add(scale_y);
    gui.add(scale_z);
    gui.add(gain);
  
  }  
};

int main() {
  MyApp app;
  app.configureAudio(48000, 512, 2, 0);
  app.start();
}
