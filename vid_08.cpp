
// Vid vis . Ars Electronica
// angle added.
// angular momentum required.
// Blood cell experiment .
// Covid is 100nm.
// Red blood cell 10,000 nm

#include "vid_include/headers_07.hpp"

using namespace al;
using namespace std;
using namespace gam;

Mesh cell_sphere, covid_sphere, spikes, back_mesh;

// Audio Parameters
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

// Audio Synth

class FM : public SynthVoice
{
public:
  // Unit generators
  gam::Pan<> mPan;
  gam::ADSR<> mAmpEnv;
  gam::ADSR<> mModEnv;
  gam::EnvFollow<> mEnvFollow;
  gam::ADSR<> mVibEnv;

  gam::Sine<> car, mod, mVib; // carrier, modulator sine oscillators

  // Additional members
  Mesh mMesh;
  float mDur;
  float mModAmt = 50;
  float mVibFrq;
  float mVibDepth;
  float mVibRise;

  void init() override
  {
    //      mAmpEnv.curve(0); // linear segments

    mAmpEnv.levels(0, 1, 1, 0);
    mModEnv.levels(0, 1, 1, 0);
    mVibEnv.levels(0, 1, 1, 0);
    //      mVibEnv.curve(0);

    // We have the mesh be a sphere
    addDisc(mMesh, 1.0, 30);
  }

  //
  void onProcess(AudioIOData &io) override
  {
    updateFromParameters();
    mVib.freq(mVibEnv());
    float carBaseFreq =
        shared_freq[covi] * carMul[covi];
    float modScale =
        shared_freq[covi] * modMul[covi];
    float gain = shared_amp[covi] * 0.1;

    while (io())
    {
      mVib.freq(mVibEnv());
      car.freq((1 + mVib() * mVibDepth) * carBaseFreq +
               mod() * mModEnv() * modScale);
      float s1 = car() * mAmpEnv() * gain;
      float s2;
      mEnvFollow(s1);
      mPan(s1, s1, s2);
      io.out(0) += s1;
      io.out(1) += s2;
    }
    if (mAmpEnv.done() && (mEnvFollow.value() < 0.001))
      free();
  }

  void onTriggerOn() override
  {
    updateFromParameters();

    float modFreq = shared_freq[covi];
    mod.freq(modFreq);

    mVibEnv.lengths()[0] = (1 - mVibRise);
    mVibEnv.lengths()[1] = mVibRise;
    mAmpEnv.reset();
    mVibEnv.reset();
    mModEnv.reset();
  }
  void onTriggerOff() override
  {
    mAmpEnv.triggerRelease();
    mModEnv.triggerRelease();
    mVibEnv.triggerRelease();
  }

  void updateFromParameters()
  {
    mModEnv.levels()[0] = idx1[covi];
    mModEnv.levels()[1] = idx2[covi];
    mModEnv.levels()[2] = idx2[covi];
    mModEnv.levels()[3] = idx3[covi];
    mod.freq(shared_freq[covi]);
    mAmpEnv.levels()[1] = 1.0;
    // mAmpEnv.levels()[2] = 1.0;
    mAmpEnv.levels()[2] = sustain[covi];

    mAmpEnv.lengths()[0] = attackTime[covi];
    mModEnv.lengths()[0] = attackTime[covi];

    mAmpEnv.lengths()[3] = releaseTime[covi];
    mModEnv.lengths()[3] = releaseTime[covi];

    // mAmpEnv.totalLength(getInternalParameterValue("dur"), 1);
    mModEnv.lengths()[1] = mAmpEnv.lengths()[1];
    mVibEnv.levels()[1] = vibRate1[covi];
    mVibEnv.levels()[2] = vibRate2[covi];
    mVibDepth = vibDepth[covi];
    mVibRise = vibRise[covi];
    mPan.pos(pan[covi]);
  }
};

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
    // col = RGB(1,1,1);
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
  Parameter carMulti{"/carMulti", "", 0.2, "", 0., 20.0};
  Parameter modMulti{"/modMulti", "", 1.01, "", 0., 20.0};
  //  Parameter gain{"gain", 0.01, 0, 1};
  Parameter frequency{"Frequency", 60, 0, 127};
  Parameter modulation{"Modulation", 60, 0, 127};
  Parameter index{"Index", 60, 0, 127};
  Parameter gain{"Gain", -90, -90, 0};
  // Boid boids[Nb];
  Mesh tails;
  Mesh box;
  Mesh mesh;
  vector<Cells> cells;
  vector<Covids> covids;
  vector<Node> node;

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
  float back_color_phase = 0.1;
  float space_scale = 150;
  // Audio Synth
  SynthGUIManager<FM> synthManager{"Vid_synth"};
  //    ParameterMIDI parameterMIDI;
  int midiNote;
  int current_nav;
  bool change_nav;
  bool follow_cell;

  bool freeze = false;
  float timer = 0;
 int frame[Nb];

  void onCreate()
  {
    current_nav = 0;
    change_nav = false;
    follow_cell = false;
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
    // Randomly distribute the starting frames
    for (int i = 0; i < Nb; i++)
    {
      // frame[i] = (int)al::rnd::uniform(999);
      frame[i] = (int)al::rnd::uniform(0);

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
    { // Spike...
      Mesh &m = spikes;
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

  void onAnimate(double dt)
  {
    imguiBeginFrame();
    synthManager.drawSynthControlPanel();
    imguiEndFrame();
    a += 0.29;
    b += 0.23;
    //    rand_vec = {scale_x.get() * al::rnd::uniform(0.5, 1.), scale_y.get() * al::rnd::uniform(0.5, 1.), scale_z.get() * al::rnd::uniform(0.5, 1.)};
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
      cells[i].col = HSV(float(i) / Nb * 0.01 + 1, 0.5, 0.6);
      // heads.color(HSV(rnd::uniform(), 1.0f, 1.0f));

      tails.vertex(cells[i].pos);
      tails.vertex(cells[i].pos - cells[i].vel.normalized(0.5)); //꼬리 길이 0.07otherwise it will lenght of velo
      tails.color(RGB(1, 0.2, 0.5));
      // Define covids position
      covids[i].pos = cells[i].pos + (float)((0.05 * data[i][frame[i]]) + cell_radius) * init_vec[i].normalized();
      covids[i].col = HSV(float(i) / Nb * 0.01 + 0.1, 0.2, 1);
      cout<< data[i][frame[i]] <<endl;
    }

    //  frequency.set(::map(data[Nb][frame], -19, 32, 0, 127));
    //   modulation.set(::map(data[Nb+1][frame], -19, 32, 0, 127));
    //   index.set(::map(data[Nb+2][frame],-19, 32, 0, 127));
    

    // Press f to follow and h to change
    if (follow_cell)
    {
      if (change_nav)
      {
        current_nav += 1;
        change_nav = !change_nav;
        if (current_nav == Nb - 1)
        {
          current_nav = 0;
        }
      }
      nav().pos(cells[current_nav].pos + Vec3f{0, 0, 10});
      nav().faceToward(nav().pos(), cells[current_nav].pos);
    }
    // cout << nav().pos() << endl;

    for (int i = 0; i < Nb; ++i)
    {
      // // Audio parameter Update. We can update synth parameters here with any data
      // // Collision detection, Number of Bonds, Binding force, Penetrate(Sustain), Tension (release, vibration)
      // // Detection algorithm required
      shared_freq[i] = data[i][frame[i]] * 50;
      // shared_amp[Nb] = sqrt(1 / (nav().pos() - covids[i].pos).mag());
      shared_amp[i] = 0.01;
      //      shared_amp[Nb] = data[i][frame];
      attackTime[i] = 0.1;
      releaseTime[i] = 2;
      sustain[i] = 3;
      idx1[i] = 0.1;
      idx2[i] = 3;
      idx3[i] = 5;
      // carMul[i] = carMulti.get();
      // modMul[i] = modMulti.get();
      carMul[i] = 1.8;
      modMul[i] = 1.5;
      vibRate1[i] = 1;
      vibRate2[i] = 1;
      vibRise[i] = 4;
      vibDepth[i] = 3;
      pan[i] = 0.5;
    }
    //     for (int i = 0; i < nodeCount; i++)
    //     {
    //       back_color_phase = back_color_phase + 0.1;
    //       Node &n = node[i];
    //       n.col = HSV(0.9, back_color_phase, 1.0);
    // //      n.position = n.position + randomVec3f(0.1);
    //     }

    for (int i = 0; i < Nb; ++i)
    {
      // Audio trigger loop
      // Event detection: Collision. Attack time.
      if (frame[i] == 0)
      {
        covi = i;
        synthManager.triggerOn(covi);
        cout << "Trigger on: " << i << endl;
      }
      //  Update audio parameters
      frame[i] += 1;
      ///////////////////
      //      cout << cells[0].angle << endl;
      if (frame[i] == 999)
      {
        covi = i;
        synthManager.triggerOff(covi);
        cout << "Trigger off: " << i << endl;
        frame[i] = 0;
      }
    }
  }

  void onDraw(Graphics &g)
  {
    g.clear(0);
    // synthManager.render(g);
    // Draw GUI
    // imguiDraw();
    // Draw
    g.lighting(10);
    g.depthTesting(true);
    g.blending(true);
    g.blendTrans();
    // g.pointSize(10); //전체 포인트 사이즈
    // g.nicest();support
    // g.stroke(10);
    // cell_sphere.reset();
    // covid_sphere.reset();
    tails.reset();
    g.rotate(a, Vec3f(0, 1, 0));
    g.rotate(b, Vec3f(1));
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
    texture.unbind();
    back_mesh.reset(); /////////////////////////////////////////////////////////////
    g.color(1);
    if (box_draw)
    {
      g.draw(box);
      g.draw(tails);
    }
  }
  void onSound(AudioIOData &io) override
  {
    synthManager.render(io); // Render audio
     while (io()) {
    //  io.out(0) = mOsc1() + mOsc3() * gain;
    //  io.out(1) = mOsc1() + mOsc2() * (1-gain);
    //   modulator.freq(mtof(modulation.get()));
    //   carrier.freq(mtof(frequency.get()) + mtof(index.get()) * modulator());
    //   float v = carrier() * dbtoa(gain.get());
    //   io.out(0) = io.out(1) = v;

    }
  }

  bool onKeyDown(const Keyboard &k)
  {
    switch (k.key())
    {
    case 'r':
      resetCells();
      break;
    case 'f':
      follow_cell = !follow_cell;
      break;
    case 'h':
      change_nav = !change_nav;
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
    gam::sampleRate(audioIO().framesPerSecond());
    auto GUIdomain = GUIDomain::enableGUI(defaultWindowDomain());
    auto &gui = GUIdomain->newGUI();
    gui.add(time);
    gui.add(separation);
    gui.add(separationStrength);
    gui.add(alignment);
    gui.add(centering);
    gui.add(huntUrge);
    gui.add(carMulti);
    gui.add(modMulti);
        gui.add(gain).add(modulation).add(index).add(frequency);

  }
};

int main() { MyApp().start(); }
