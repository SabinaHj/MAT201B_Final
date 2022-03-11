
// Vid vis . Ars Electronica
// angle added. 
// angular momentum required.
// Blood cell experiment . 
//Covid is 100nm. 
//Red blood cell 10,000 nm

#include <cmath>
#include "al/app/al_App.hpp"
#include "al/graphics/al_Shapes.hpp"
#include "al/math/al_Functions.hpp"
#include "al/math/al_Random.hpp"
#include "al/app/al_GUIDomain.hpp"
#include "Gamma/Noise.h"
#include "Gamma/Delay.h"
#include "Gamma/Oscillator.h"
#include "Gamma/Filter.h"
#include "Gamma/Envelope.h"
#include "Gamma/Effects.h"
#include "Gamma/Analysis.h"
#include "al/io/al_CSVReader.hpp"

#include <fstream>
#include <vector>
using namespace al;
using namespace std;
using namespace gam;

Mesh cell_sphere, covid_sphere, spikes;
Pan<> mPan;
Sine<> mOsc1;
Sine<> mOsc2;
Sine<> mOsc3;
Env<3> mAmpEnv;
EnvFollow<> mEnvFollow;
NoisePink<> noise;
Comb<float, ipl::Switchable> comb;
double a = 0;
double b = 0;
LFO<> mod;
OnePole<> onePole;
Vec3f rand_vec;

struct Cells
{
  Vec3f pos, vel, acc, axis, vib, ang_moment;
  Color col;
  float angle;
  void update(float dt)
  {
    pos += vel * dt;
    // float dt = dt_ms;
  }
  void draw(Graphics &g)
  {
    g.pushMatrix();
    g.translate(pos);
    g.rotate(angle, vel);
    g.scale(vib);
    g.color(col);
    g.draw(cell_sphere);
    g.popMatrix();
  }
};

struct Covids
{
  Vec3f pos, vel, acc, axis, vib, ang_moment;
  Color col;
  float angle;
  void update(float dt)
  {
    pos += vel * dt;
    // float dt = dt_ms;
  }
  void draw(Graphics &g)
  {
    g.pushMatrix();
    g.translate(pos);
//    g.rotate(angle, vel);
//    g.scale(vib);
    g.color(col);
    g.draw(covid_sphere);
    g.popMatrix();
  }
};

// float pushRadius = 0.05;
// float matchRadius = 0.125;

struct MyApp : public App
{
  Parameter time{"/time", "", 0.01, "", 0.001, 0.1};
  Parameter separation{"/separation", "", 0.02, "", 0.01, 1.0};
  Parameter separationStrength{"/separationStrength", "", 0.05, "", 0.01, 1.0};
  Parameter alignment{"/alignment", "", 0.125, "", 0.001, 0.3};
  Parameter centering{"/centering", "", 0.05, "", 0.01, 0.2};
  Parameter huntUrge{"/huntUrge", "", 0.2, "", 0.01, 1};
  Parameter scale_x{"x", 1, 0, 2};
  Parameter scale_y{"y", 1, 0, 2};
  Parameter scale_z{"z", 1, 0, 2};
  Parameter gain{"gain", 0.01, 0, 1};

  static const int Nb = 10; // Number of boids
  // Boid boids[Nb];
  Mesh tails;
  Mesh box;
  Mesh mesh;
  vector<Cells> cells;
  vector<Covids> covids;

  float force;
  bool box_draw = false;
  Vec3f unitVector;
  Vec3f gravity;
  int s;
  float cell_angle;
// outside is 422 inside is 1000 each
  float data[422][1000];

  void onCreate()
  {
    CSVReader reader;
    reader.addType(CSVReader::REAL);
    reader.readFile("data/Y_voltage_force_flatten_transpose.csv");
    std::vector<double> column0 = reader.getColumn(0);
    for (int i = 0; i < 422; i++) {
        for (int j = 0; j < 1000; j++) {
            data[i][j] = column0[j+i*1000];
        }
    }

    cell_angle = al::rnd::uniform();
    nav().pos(0.5, 0.7, 20);
    nav().faceToward(Vec3d(0, 0, 0), Vec3d(0, 1, 0));
    addSphere(cell_sphere, 2, 50, 50);
    addSphere(covid_sphere, 1, 50, 50);

 //   blub spikes
    addSurface(spikes,
               33,  // number of points along x
               33   // number of points along y
    );
    {
      Mesh &m = spikes;
      for (int i = 0; i < m.vertices().size(); ++i) {
        Mesh::Vertex &v = m.vertices()[i];
        float r = ::hypot(v.x, v.y);
        v.z = ::exp(-8 * r * r);
      }
    }

    cell_sphere.generateNormals();
    covid_sphere.generateNormals();
    cells.resize(Nb);
    covids.resize(Nb);

    // box.primitive(Mesh::LINE_LOOP);
    box.primitive(Mesh::LINES);
    s = 10;
    box.vertex(s, s, s);
    box.vertex(-s, s, s);
    box.vertex(s, s, s);
    box.vertex(s, -s, s);
    box.vertex(s, s, s);
    box.vertex(s, s, -s);

    box.vertex(s, -s, -s);
    box.vertex(-s, -s, -s);
    box.vertex(s, -s, -s);
    box.vertex(s, s, -s);
    box.vertex(s, -s, -s);
    box.vertex(s, -s, s);

    box.vertex(-s, s, -s);
    box.vertex(s, s, -s);
    box.vertex(-s, s, -s);
    box.vertex(-s, -s, -s);
    box.vertex(-s, s, -s);
    box.vertex(-s, s, s);

    box.vertex(-s, -s, s);
    box.vertex(s, -s, s);
    box.vertex(-s, -s, s);
    box.vertex(-s, s, s);
    box.vertex(-s, -s, s);
    box.vertex(-s, -s, -s);

    resetCells(); // reset
  }



  // Randomize boid positions/velocities uniformly inside unit disc
  void resetCells()
  {
    for (auto &c : cells)
    {
      c.pos = al::rnd::ball<Vec3f>()*5;
      c.vel = al::rnd::ball<Vec3f>();
    }
  }

  bool freeze = false;
  float timer = 0;
  int frame = 0;

  void onAnimate(double dt)
  {
   frame += 1;
   a += 0.29;
   b += 0.23;
   rand_vec = {scale_x.get() * al::rnd::uniform(0.5, 1.), scale_y.get() * al::rnd::uniform(0.5, 1.), scale_z.get() * al::rnd::uniform(0.5, 1.)};
    if (freeze)
      return;

    dt = time;
    timer +=dt;
    //     cout << dt << endl;
    // Compute boid-boid interactions
    // vector<Vec3f> &position(mesh.vertices());
    for (int i = 0; i < Nb; ++i)
    {
      for (int j = i + 1; j < Nb; ++j)
      {
        // printf("checking cells %d and %d\n", i,j);

        auto ds = cells[j].pos - cells[i].pos;
        auto dist = ds.mag();

        // dist = position[j] - position[i];

        // Collision avoidance
        // float pushRadius = 0.05;
        // float pushStrength = 1;
        float push = exp(-al::pow2(dist / separation)) * separationStrength;

        auto pushVector = ds.normalized() * push;
        cells[i].vel += pushVector;
        cells[j].vel -= pushVector;

        // Velocity matching
        // float matchRadius = 0.125;
        float nearness = exp(-al::pow2(dist / alignment));
        Vec3f veli = cells[i].vel;
        Vec3f velj = cells[j].vel;

        // Take a weighted average of velocities according to nearness
        cells[i].vel = veli * (1 - 0.5 * nearness) + velj * (0.5 * nearness);
        cells[j].vel = velj * (1 - 0.5 * nearness) + veli * (0.5 * nearness);
        //cells[i].ang_moment = cells[j].vel.mag()*10; 
        cells[i].angle += cells[j].vel.mag()*100;
        cells[i].vib = Vec3f{0.7, 0.7, 0.7} + cells[j].vel* 0.1;
      }
    }

    // TODO: Flock centering
    // float pushRadius = 0.05;
    // float matchRadius = 0.125;

    for (int j = 0; j < Nb; j++)
    {
      // find the center of the flock
      const Vec3f &point = cells[j].pos;

      int count = 0;
      Vec3f sum(0, 0, 0);
      cells[j].acc = {0, 0, 0};
      for (int i = 0; i < Nb; i++)
      {
        if ((point - cells[i].pos).mag() < separation)
        {
          sum += cells[i].pos; // find the center
          count++;
          // find the bounds of these points
          // cout << count << endl;
        }
      }
      if (count > 0)
      {
        Vec3f center = sum / count; // what if count == 0?
        auto dsc = cells[j].pos - center;
        if (dsc.magSqr() > 0)
        {
          force = centering / dsc.magSqr();
          if (force > 1)
          {
            force = 1;
          }
          unitVector = dsc.normalized();
          gravity = unitVector * force;
          cells[j].acc += gravity;
          cells[j].vel += cells[j].acc;
        }
      }
    }

    // for (int i = 0; i < Nb; i++){
    //   // pull i toward the center
    //   // assignment 2 pull toward == gravity
    // }

    // Update boid independent behaviors
    for (auto &b : cells)
    {
      // Random "hunting" motion
      // float huntUrge = 0.5;
      auto hunt = al::rnd::ball<Vec3f>();
      // Use cubed distribution to make small jumps more frequent
      hunt *= hunt.magSqr();
      b.vel += hunt * huntUrge;

      // wraparound
      int r = s;
//      int d = r * 2;

      if (b.pos.x > r)
      {
        b.vel.x = -b.vel.x;
      }
      if (b.pos.x < -r)
      {
        b.vel.x = -b.vel.x;
      }
      if (b.pos.y > r)
      {
        b.vel.y = -b.vel.y;
      }
      if (b.pos.y < -r)
      {
        b.vel.y = -b.vel.y;
      }
      if (b.pos.z > r)
      {
        b.vel.z = -b.vel.z;
      }
      if (b.pos.z < -r)
      {
        b.vel.z = -b.vel.z;
      }
    }

    // Generate meshes

    // heads.reset();//reset clears all
    //    heads.primitive(Mesh::POINTS);
    // heads.primitive(Mesh::POINTS);
    //    cell_sphere.reset();
    tails.reset();
    tails.primitive(Mesh::LINES);

    for (int i = 0; i < Nb; ++i)
    {
      cells[i].update(dt);
      // heads.vertex(cells[i].pos);
      // heads.color(HSV(float(i) / Nb * 0.3 + 0.3, 0.7));
      cells[i].col = HSV(float(i) / Nb * 0.01 + 1, 0.2, 0.6);
      // heads.color(HSV(rnd::uniform(), 1.0f, 1.0f));

      tails.vertex(cells[i].pos);
      tails.vertex(cells[i].pos - cells[i].vel.normalized(0.5)); //꼬리 길이 0.07otherwise it will lenght of velo
      tails.color(RGB(1, 0.2, 0.5));
      // Define covids position
      covids[i].pos = cells[i].pos + (0.5* data[i][frame]);
      covids[i].col = HSV(float(i) / Nb * 0.01 + 0.1, 0.2, 0.6);
    }
//      cout << cells[0].angle << endl;    
    if (frame == 999){
      frame = 0;
    }
  }

  void onDraw(Graphics &g)
  {
    g.clear(0);
    g.lighting(true);
    g.depthTesting(true);
    // g.pointSize(10); //전체 포인트 사이즈
    // g.nicest();support
    // g.stroke(10);
    g.rotate(a, Vec3f(0, 1, 0));
    g.rotate(b, Vec3f(1));
    // g.scale(rand_vec);
    g.meshColor();

    // Cells
    for (int i = 0; i < Nb; i++)
    {
      g.pushMatrix();
      // texture.bind();
      for (auto c : cells)
        c.draw(g);
        g.scale(rand_vec);
      // texture.unbind();
      g.popMatrix();
    }

    // Covids
    for (int i = 0; i < Nb; i++)
    {
      g.pushMatrix();
      // texture.bind();
      for (auto v : covids)
        v.draw(g);
        g.scale(rand_vec);
      // texture.unbind();
      g.popMatrix();
    }

    // g.draw(mesh);

    // g.stroke(1);
    g.color(1);
    if (box_draw)
    {
      g.draw(box);
      g.draw(tails);
    }
  }

  void onSound(AudioIOData &io) override {
    mOsc1.freq(rand_vec[0]*80+50);
    mOsc2.freq(rand_vec[1]*500+10);
    mOsc3.freq(rand_vec[2]*1000);

    // mOsc1.freq(scale_x.get()*100+50);
    // mOsc2.freq(scale_y.get()*1000+50);
    // mOsc3.freq(scale_z.get()*3000+50);


//		mod.period(scale_x.get());
//		comb.delay(mod.triU() * scale_y.get() /1000 + scale_z.get()/10);
    while (io()) {
//      io.out(0) = mOsc1() + mOsc3() * gain;
 //     io.out(1) = mOsc1() + mOsc2() * (1-gain);

    }
  }

  bool onKeyDown(const Keyboard &k)
  {
    switch (k.key())
    {
    case 'r':
      resetCells();
      break;

    case ']':
      box_draw = !box_draw;
      break;

    case ' ':
      freeze = !freeze;
      break;
    }

    return true;
  }
  void onInit() override
  {
    auto GUIdomain = GUIDomain::enableGUI(defaultWindowDomain());
    auto &gui = GUIdomain->newGUI();
    gui.add(time);
    gui.add(separation);
    gui.add(separationStrength);
    gui.add(alignment);
    gui.add(centering);
    gui.add(huntUrge);
    gui.add(scale_x);
    gui.add(scale_y);
    gui.add(scale_z);
    gui.add(gain);
  
  }
};

int main() { MyApp().start(); }

/*
/ TODO: Flock centering
    for (int j = 0; j < Nb; j++) {
      // find the center of the flock
      const Vec3f& point = boids[j].pos;

      int count = 0;
      Vec3f sum(0, 0, 0);
      for (int i = 0; i < Nb; i++) {
        if ((point - boids[i].pos).mag() < radius) {
          sum += boids[i].pos;  // find the center
          count++;
          // find the bounds of these points
        }
      }
      Vec3f center = sum / count;  // what if count == 0?

      // pull i toward the center
      // assignment 2 pull toward == gravity
    }
    */
