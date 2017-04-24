using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Runtime.InteropServices;

namespace PBF_GPU_SLOW_2D {
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
            Gizmos.DrawWireCube(range / 2, range);
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
            if (!random_start) {
                for (int i = 1; i < range.x / 2; i++) {
                    for (int j = (int)(range.y / 4); j < range.y - 1; j++) {
                        particles.Add(new Particle(new Vector2(i, j), 1f));
                    }
                }
            } else {
                for (int i = 0; i < 30000; i++) {
                    particles.Add(new Particle(new Vector2(Random.value * (float)range.x/2, Random.value * (float)range.y), 1));
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
        public Vector2 oldPos;
        public Vector2 newPos;
        public Vector2 velocity;
        public Vector2 force;
        public Vector2 deltaP;
        public float mass;
        public float lambda;
        public float pConstraint;

        public Particle(Vector2 pos, float mass) {
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
