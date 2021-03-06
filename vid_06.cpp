
// Vid vis . Ars Electronica
// angle added.
// angular momentum required.
// Blood cell experiment .
// Covid is 100nm.
// Red blood cell 10,000 nm

#include "vid_include/headers_06.hpp"

using namespace al;
using namespace std;
using namespace gam;

Mesh cell_sphere, covid_sphere, spikes, back_mesh;
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
float halfSize = 1.;

LFO<> mod;
OnePole<> onePole;
Vec3f rand_vec;

Vec3f randomVec3f(float scale)
{
  return Vec3f(al::rnd::uniformS(), al::rnd::uniformS(), al::rnd::uniformS()) * scale;
}

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
  void draw(Graphics &g, Mesh &m)
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
  void draw(Graphics &g, Mesh &m)
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

struct Node
{
  Vec3f position;
  Color col;
  Node()
  {
    position = randomVec3f(150);
    //   vel =  Vec3f( 0, (rnd::uniform()-0.5) , 0);
    //col = RGB(1,1,1);
  }
  void draw(Graphics &g, Mesh &m)
  {

    m.vertex(position);
    m.color(col);
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
//  Parameter gain{"gain", 0.01, 0, 1};

  Parameter frequency{"Frequency", 100, 0, 1000};
  Parameter modulation{"Modulation", 60, 0, 127};
  Parameter index{"Index", 60, 0, 127};
  Parameter gain{"Gain", -90, -90, 0};


  static const int Nb = 10; // Number of boids
  // Boid boids[Nb];
  Mesh tails;
  Mesh box;
  Mesh mesh;
  vector<Cells> cells;
  vector<Covids> covids;
  ShaderProgram shader;
  float vid_radius = 1;
  float cell_radius = 3;
  float force;
  bool box_draw = false;
  Vec3f unitVector;
  Vec3f gravity;
  Vec3f init_vec[Nb];
  int s, boundary;
  float sounda[422];
  float cell_angle;
  // outside is 422 inside is 1000 each
  float data[422][1000];
  Texture texture;
  vector<Node> node;
  float back_color_phase = 0.1;

  gam::Sine<> carrier;
  gam::Sine<> modulator;

  void onCreate()
  {
    boundary = 60;
    node.resize(nodeCount);

    CSVReader reader;
    reader.addType(CSVReader::REAL);
    reader.readFile("data/Y_voltage_force_flatten_transpose.csv");
    std::vector<double> column0 = reader.getColumn(0);
    for (int i = 0; i < 422; i++)
    {
      for (int j = 0; j < 1000; j++)
      {
        data[i][j] = column0[j + i * 1000];
      }
    }

    cell_angle = al::rnd::uniform();
    nav().pos(0.5, 0.7, 20);
    nav().faceToward(Vec3d(0, 0, 0), Vec3d(0, 1, 0));
    addSphere(cell_sphere, cell_radius, 30, 30);
    addSphere(covid_sphere, vid_radius, 30, 30);
    back_mesh.primitive(Mesh::POINTS);
    shader.compile(vertex, fragment, geometry);

    for (int i = 0; i < Nb; ++i)
    {
      init_vec[i] = randomVec3f(1.0);
    }
    //   blub spikes
    addSurface(spikes,
               33, // number of points along x
               33  // number of points along y
    );
    {
      Mesh &m = spikes;
      for (int i = 0; i < m.vertices().size(); ++i)
      {
        Mesh::Vertex &v = m.vertices()[i];
        float r = ::hypot(v.x, v.y);
        v.z = ::exp(-8 * r * r);
      }
    }
    for (int i = 0; i < nodeCount; i++)
    {
      Node &n = node[i];
      n.col = HSV(0.9, back_color_phase, 1.0);
      //      n.position = Vec3f((al::rnd::uniform() - 0.5) * space_scale, (al::rnd::uniform() - 0.5) * space_scale, (al::rnd::uniform() - 0.5) * space_scale);
    }

    cell_sphere.generateNormals();
    covid_sphere.generateNormals();
    cells.resize(Nb);
    covids.resize(Nb);

    texture.create2D(300, 300, Texture::R8, Texture::RED, Texture::SHORT);
    int Nx = texture.width();
    int Ny = texture.height();
    std::vector<short> alpha;
    alpha.resize(Nx * Ny);
    for (int j = 0; j < Ny; ++j)
    {
      float y = float(j) / (Ny - 1) * 2 - 1;
      for (int i = 0; i < Nx; ++i)
      {
        float x = float(i) / (Nx - 1) * 2 - 1;
        float m = exp(-13 * (x * x + y * y));
        m *= pow(2, 15) - 1; // scale by the largest positive short int
        alpha[j * Nx + i] = m;
      }
    }
    texture.submit(&alpha[0]);

    // box.primitive(Mesh::LINE_LOOP);
    box.primitive(Mesh::LINES);
    s = boundary;
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
      c.pos = al::rnd::ball<Vec3f>() * 20;
      c.vel = al::rnd::ball<Vec3f>();
    }
  }

  bool freeze = false;
  float timer = 0;
  int frame = 0;

  void onAnimate(double dt)
  {
    frame += 1;
    a += 0.029;
    b += 0.023;
    rand_vec = {scale_x.get() * al::rnd::uniform(0.5, 1.), scale_y.get() * al::rnd::uniform(0.5, 1.), scale_z.get() * al::rnd::uniform(0.5, 1.)};
    if (freeze)
      return;

    dt = time;
    timer += dt;
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
        // cells[i].ang_moment = cells[j].vel.mag()*10;
        cells[i].angle += cells[j].vel.mag() * 100;
        cells[i].vib = Vec3f{0.7, 0.7, 0.7} + cells[j].vel * 0.01;
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
    //    tails.reset();
    tails.primitive(Mesh::LINES);

    for (int i = 0; i < Nb; ++i)
    {
      cells[i].update(dt);
      // heads.vertex(cells[i].pos);
      // heads.color(HSV(float(i) / Nb * 0.3 + 0.3, 0.7));
      cells[i].col = HSV(float(i) / Nb * 0.01 + 1, 0.2, 0.6);
      // heads.color(HSV(rnd::uniform(), 1.0f, 1.0f));

      tails.vertex(cells[i].pos);
      tails.vertex(cells[i].pos - cells[i].vel.normalized(0.5)); //?????? ?????? 0.07otherwise it will lenght of velo
      tails.color(RGB(1, 0.2, 0.5));
      // Define covids position
      covids[i].pos = cells[i].pos + (float)((0.05 * data[i][frame]) + cell_radius) * init_vec[i].normalized();
      covids[i].col = HSV(float(i) / Nb * 0.01 + 0.1, 0.2, 0.3);
    }
    //      cout << cells[0].angle << endl;
    if (frame == 999)
    {
      frame = 0;
    }
  }

  void onDraw(Graphics &g)
  {
    g.clear(0);

    g.lighting(true);
    g.depthTesting(true);
    g.blending(true);
    g.blendTrans();
    // g.pointSize(10); //?????? ????????? ?????????
    // g.nicest();support
    // g.stroke(10);
    // cell_sphere.reset();
    // covid_sphere.reset();
    tails.reset();
    g.rotate(a, Vec3f(0, 1, 0));
    g.rotate(b, Vec3f(1));
    // g.scale(rand_vec);
    g.meshColor();

    // Cells
    //    cell_sphere.reset();
    for (int i = 0; i < Nb; i++)
    {
      g.pushMatrix();
      for (auto c : cells)
      {
        c.draw(g, cell_sphere);
      }
      g.shader(shader);
      g.popMatrix();
    }

    // Covids
    //    covid_sphere.reset();
    for (int i = 0; i < Nb; i++)
    {
      g.pushMatrix();
      for (auto v : covids)
      {
        v.draw(g, covid_sphere);
      }
      g.popMatrix();
    }

    back_mesh.reset(); /////////////////////////////////////////////////////////////
    for (int i = 0; i < nodeCount; i++)
    {
      texture.bind();
      for (auto n : node)
        n.draw(g, back_mesh);
    }
    g.meshColor();
    g.shader(shader);
    g.shader().uniform("halfSize", halfSize);
    g.draw(back_mesh);
    texture.unbind();
    back_mesh.reset(); /////////////////////////////////////////////////////////////

    //    g.draw(cell_sphere);
    //    g.draw(covid_sphere);
    // g.draw(mesh);
    // g.stroke(1);
    g.color(1);
    if (box_draw)
    {
      g.draw(box);
      g.draw(tails);
    }
  }
  int i = 0;
  int k = 0;
  float asdf =0;
  void onSound(AudioIOData &io) override
  {

    while (io())
    {
//      modulator.freq(mtof(modulation.get()));
      for (int vid = 0; vid<1; vid++){
        carrier.freq(frequency.get() + data[vid][frame]);
        float v = carrier() * dbtoa(gain.get());
        sounda[vid] = v * 1000;
        // cout<< "    " << v << "         " <<sounda << endl;
      }
      for (int vid = 0; vid<10; vid++){
        asdf += sounda[vid];
      }
      io.out(0) = io.out(1) = asdf;

      cout << frame << "    " << asdf << endl;

      // cout << i << "  " << k << "    " <<  sounda << endl;
      
      // k++;
      // if (k > 200)
      // {
      //   i++;
      //   k = 0;
      // }

      // if (i == 999)
      // {
      //   i = 0;
      // }
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
    case '.':
      boundary++;
      cout << "boundary: " << boundary << endl;
      break;
    case ',':
      boundary--;
      cout << "boundary: " << boundary << endl;
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
    // gui.add(scale_x);
    // gui.add(scale_y);
    // gui.add(scale_z);
    gui.add(gain);
    gui.add(frequency);
    gui.add(modulation);
    gui.add(index);
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
