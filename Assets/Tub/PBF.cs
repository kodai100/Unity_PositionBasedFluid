using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Runtime.InteropServices;

namespace PBF_GPU_SLOW_2D_TUB {
    public class PBF : MonoBehaviour {

        List<Particle> particles = new List<Particle>();
        Particle[] particle_array;
        [SerializeField]int maxParticleNum;

        public Vector2 GRAVITY = new Vector2(0f, -9.8f);
        public int PRESSURE_ITERATIONS = 2;
        public float dt = 0.005f;
        public Vector2 range = new Vector2(30, 30);
        public bool random_start = false;

        #region GPU
        public ComputeShader PBF_CS;
        ComputeBuffer particleBufferRead;
        ComputeBuffer particleBufferWrite;
        static int SIMULATION_BLOCK_SIZE = 256;
        int threadGroupSize;
        #endregion GPU


        public float H = 1.2f;
        float KPOLY;
        float SPIKY;
        float VISC;
        public float REST_DENSITY = 1f;

        static float EPSILON_LAMBDA = 150f;
        public float C = 0.01f;  // Typically 0.01

        static float EPSILON_VORTICITY = 10f;
        static float K = 0.001f;
        static float deltaQMag;
        static float wQH;

        #region MonoBehaviour

        void Start() {
            KPOLY = (float)(315f / (64f * Mathf.PI * Mathf.Pow(H, 9)));
            SPIKY = (float)(45f / (Mathf.PI * Mathf.Pow(H, 6)));
            VISC = (float)(15f / (2 * Mathf.PI * (H * H * H)));
            deltaQMag = 0.3f * H;
            wQH = KPOLY * (H * H - deltaQMag * deltaQMag) * (H * H - deltaQMag * deltaQMag) * (H * H - deltaQMag * deltaQMag);

            CreateWater();
            InitializeComputeBuffer();
        }

        void Update() {

            PositionBasedFluid();

        }

        void OnDestroy() {
            if(particleBufferRead != null) {
                particleBufferRead.Release();
                particleBufferRead = null;
            }
            if(particleBufferWrite != null) {
                particleBufferWrite.Release();
                particleBufferWrite = null;
            }
        }

        void OnDrawGizmos() {
            // Gizmos.DrawWireCube(range / 2, range);
        }

        #endregion MonoBehaviour

        #region GPUProcess

        void InitializeComputeBuffer() {
            particleBufferRead = new ComputeBuffer(maxParticleNum, Marshal.SizeOf(typeof(Particle)));
            particleBufferWrite = new ComputeBuffer(maxParticleNum, Marshal.SizeOf(typeof(Particle)));
            particleBufferRead.SetData(particle_array);
            particleBufferWrite.SetData(particle_array);
            threadGroupSize = Mathf.CeilToInt(maxParticleNum / SIMULATION_BLOCK_SIZE) + 1;
        }

        void PositionBasedFluid() {

            // Debug.Log(maxParticleNum);

            // Set Variables
            PBF_CS.SetInt("_NumParticles", maxParticleNum);
            PBF_CS.SetVector("_Gravity", GRAVITY);
            PBF_CS.SetFloat("_DT", dt);
            PBF_CS.SetFloat("_H", H);
            PBF_CS.SetFloat("_KPoly", KPOLY);
            PBF_CS.SetFloat("_Spiky", SPIKY);
            PBF_CS.SetFloat("_Visc", VISC);
            PBF_CS.SetFloat("_RestDensity", REST_DENSITY);
            PBF_CS.SetFloat("_EpsilonLambda", EPSILON_LAMBDA);
            PBF_CS.SetFloat("_C", C);
            PBF_CS.SetFloat("_EpsilonVorticity", EPSILON_VORTICITY);
            PBF_CS.SetFloat("_K", K);
            PBF_CS.SetFloat("_DeltaQMag", deltaQMag);
            PBF_CS.SetFloat("_WQH", wQH);
            PBF_CS.SetVector("_Range", range);

            

            int kernel;
            // Copy Buffer Every Frame
            kernel = PBF_CS.FindKernel("CopyBuffer");
            PBF_CS.SetBuffer(kernel, "_ParticlesBufferRead", particleBufferRead);
            PBF_CS.SetBuffer(kernel, "_ParticlesBufferWrite", particleBufferWrite);
            PBF_CS.Dispatch(kernel, threadGroupSize, 1, 1);

            // Process1. Update
            kernel = PBF_CS.FindKernel("Update");
            PBF_CS.SetBuffer(kernel, "_ParticlesBufferRead", particleBufferRead);
            PBF_CS.SetBuffer(kernel, "_ParticlesBufferWrite", particleBufferWrite);
            PBF_CS.Dispatch(kernel, threadGroupSize, 1, 1);

            // Process2 Lambda Iteration
            for (int i = 0; i < PRESSURE_ITERATIONS; i++) {

                // CalcLambda
                kernel = PBF_CS.FindKernel("CalcLambda");
                PBF_CS.SetBuffer(kernel, "_ParticlesBufferRead", particleBufferRead);
                PBF_CS.SetBuffer(kernel, "_ParticlesBufferWrite", particleBufferWrite);
                PBF_CS.Dispatch(kernel, threadGroupSize, 1, 1);

                // CalcDeltaP
                kernel = PBF_CS.FindKernel("CalcDeltaP");
                PBF_CS.SetBuffer(kernel, "_ParticlesBufferRead", particleBufferRead);
                PBF_CS.SetBuffer(kernel, "_ParticlesBufferWrite", particleBufferWrite);
                PBF_CS.Dispatch(kernel, threadGroupSize, 1, 1);

                // ApplyDeltaP
                kernel = PBF_CS.FindKernel("ApplyDeltaP");
                PBF_CS.SetBuffer(kernel, "_ParticlesBufferRead", particleBufferRead);
                PBF_CS.SetBuffer(kernel, "_ParticlesBufferWrite", particleBufferWrite);
                PBF_CS.Dispatch(kernel, threadGroupSize, 1, 1);
            }

            // Update Velocity
            kernel = PBF_CS.FindKernel("UpdateVelocity");
            PBF_CS.SetBuffer(kernel, "_ParticlesBufferRead", particleBufferRead);
            PBF_CS.SetBuffer(kernel, "_ParticlesBufferWrite", particleBufferWrite);
            PBF_CS.Dispatch(kernel, threadGroupSize, 1, 1);

            SwapBuffer(ref particleBufferRead, ref particleBufferWrite);
        }

        void CreateWater() {
            
            // 液体
            for (int i = 0; i < 10000; i++) {
                particles.Add(new Particle(false, new Vector2((float)range.x / 4 + Random.value * (float)range.x/2, range.y / 4 + Random.value * (float)range.y*2), 1));
            }

            // 壁
            int num_layer = 5;
            float radius = 0.2f;
            float diameter = radius * 2;
            Vector3 pos = new Vector3(0, 0, 0);

            int num_height = (int)(range.y / diameter);
            int num_width = (int)(range.x / diameter);

            // 下の壁
            for (int l = 0; l < num_layer; l++) {
                for (int i = 0; i < num_width; i++) {
                    particles.Add(new Particle(true, new Vector2(i * diameter + radius, -diameter * l - radius), 1));
                }
            }

            // 左の壁
            for (int l = 0; l < num_layer; l++) {
                for (int i = 0; i < num_height; i++) {
                    particles.Add(new Particle(true, new Vector2(diameter * l + radius, i * diameter + radius), 1));
                }
            }

            // 右の壁
            for (int l = 0; l < num_layer; l++) {
                for (int i = 0; i < num_height; i++) {
                    Vector2 pivot = new Vector2((num_width - num_layer) * diameter, pos.y);
                    particles.Add(new Particle(true, pivot + new Vector2(diameter * l + radius, i * diameter + radius), 1));
                }
            }

            particle_array = particles.ToArray();
            maxParticleNum = particle_array.Length;
        }

        void SwapBuffer(ref ComputeBuffer src, ref ComputeBuffer dst) {
            ComputeBuffer tmp = src;
            src = dst;
            dst = tmp;
        }

        #endregion GPUProcess

        public ComputeBuffer GetBuffer() {
            return particleBufferRead;
        }

        public int GetMaxParticleNum() {
            return maxParticleNum;
        }
    }


    struct Particle {
        public bool wall;
        public Vector2 oldPos;
        public Vector2 newPos;
        public Vector2 velocity;
        public Vector2 force;
        public Vector2 deltaP;
        public float mass;
        public float lambda;
        public float pConstraint;

        public Particle(bool wall, Vector2 pos, float mass) {
            this.wall = wall;
            this.oldPos = pos;
            this.mass = mass;
            this.lambda = 0;
            this.pConstraint = 0;
            this.newPos = new Vector2(0f, 0f);
            this.velocity = new Vector2(0f, 0f);
            this.force = new Vector2(0f, 0f);
            this.deltaP = new Vector2(0f, 0f);
        }
    }

}
