
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

#include <fstream>
#include <vector>
using namespace al;
using namespace std;

Mesh cell_sphere, spikes;

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

  static const int Nb = 50; // Number of boids
  // Boid boids[Nb];
  Mesh tails;
  Mesh box;
  Mesh mesh;
  vector<Cells> cells;

  float force;
  bool box_draw = false;
  Vec3f unitVector;
  Vec3f gravity;
  int s;
  float cell_angle;
  void onCreate()
  {
    cell_angle = rnd::uniform();
    nav().pos(0.5, 0.7, 20);
    nav().faceToward(Vec3d(0, 0, 0), Vec3d(0, 1, 0));
    addSphere(cell_sphere, 0.2, 50, 50);

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
    cells.resize(Nb);

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
      c.pos = rnd::ball<Vec3f>()*5;
      c.vel = rnd::ball<Vec3f>();
    }
  }

  bool freeze = false;
  float timer = 0;
  void onAnimate(double dt)
  {
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
      auto hunt = rnd::ball<Vec3f>();
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
    }
//      cout << cells[0].angle << endl;    
  }

  void onDraw(Graphics &g)
  {
    g.clear(0);
    g.lighting(true);
    g.depthTesting(true);
    // g.pointSize(10); //전체 포인트 사이즈
    // g.nicest();support
    // g.stroke(10);
    g.meshColor();

    // Finger
    for (int i = 0; i < Nb; i++)
    {
      g.pushMatrix();
      // texture.bind();
      for (auto c : cells)
        c.draw(g);
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
