using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Runtime.InteropServices;

namespace PBF_GPU_FAST_2D {
    public class PBF : MonoBehaviour {

        List<Particle> particles = new List<Particle>();
        Particle[] particle_array;

        [SerializeField]int maxParticleNum;
        public enum Mode {
            NUM_8K, NUM_16K, NUM_32K, NUM_65K, NUM_130K, NUM_260K
        }
        public Mode num;

        public bool showGrid = false;

        public Vector2 GRAVITY = new Vector2(0f, -9.8f);
        public int PRESSURE_ITERATIONS = 4;
        public float dt = 0.005f;
        public Vector2 GridDim = new Vector2(64,64);
        public Vector2 range = new Vector2(256, 256);

        #region GPU
        public ComputeShader PBF_CS;
        ComputeBuffer particleBufferRead;
        ComputeBuffer particleBufferWrite;
        ComputeBuffer sortedParticleBuffer;
        ComputeBuffer gridBuffer;
        ComputeBuffer gridPingPongBuffer;
        ComputeBuffer gridIndicesBuffer;
        static int SIMULATION_BLOCK_SIZE = 32;
        float gridH;
        int threadGroupSize;

        public ComputeShader BitonicCS;
        static uint BITONIC_BLOCK_SIZE = 512;
        static uint TRANSPOSE_BLOCK_SIZE = 16;
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
            gridH = range.x / GridDim.x;

            switch (num) {
                case Mode.NUM_8K:
                    maxParticleNum = 8192;
                    break;
                case Mode.NUM_16K:
                    maxParticleNum = 16384;
                    break;
                case Mode.NUM_32K:
                    maxParticleNum = 32768;
                    break;
                case Mode.NUM_65K:
                    maxParticleNum = 65536;
                    break;
                case Mode.NUM_130K:
                    maxParticleNum = 131072;
                    break;
                case Mode.NUM_260K:
                    maxParticleNum = 262144;
                    break;
                default:
                    maxParticleNum = 8192;
                    break;
            }

            InitializeComputeBuffer();
            CreateWater();
            
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
            if (gridBuffer != null) {
                gridBuffer.Release();
                gridBuffer = null;
            }
            if (gridIndicesBuffer != null) {
                gridIndicesBuffer.Release();
                gridIndicesBuffer = null;
            }
            if (gridPingPongBuffer != null) {
                gridPingPongBuffer.Release();
                gridPingPongBuffer = null;
            }
            if (sortedParticleBuffer != null) {
                sortedParticleBuffer.Release();
                sortedParticleBuffer = null;
            }
        }

        void OnDrawGizmos() {

            Gizmos.DrawWireCube(range / 2, range);

            if (showGrid) {
                Gizmos.color = Color.blue;
                for (int i = 1; i < GridDim.y; i++) {
                    Gizmos.DrawLine(new Vector3(0, gridH * i, 0), new Vector3(range.x, gridH * i, 0));
                }

                for (int i = 1; i < GridDim.x; i++) {
                    Gizmos.DrawLine(new Vector3(gridH * i, 0, 0), new Vector3(gridH * i, range.y, 0));
                }
            }
        }

        #endregion MonoBehaviour

        #region GPUProcess

        void InitializeComputeBuffer() {
            particleBufferRead = new ComputeBuffer(maxParticleNum, Marshal.SizeOf(typeof(Particle)));
            particleBufferWrite = new ComputeBuffer(maxParticleNum, Marshal.SizeOf(typeof(Particle)));
            sortedParticleBuffer = new ComputeBuffer(maxParticleNum, Marshal.SizeOf(typeof(Particle)));
            gridBuffer = new ComputeBuffer(maxParticleNum, Marshal.SizeOf(typeof(Uint2)));
            gridPingPongBuffer = new ComputeBuffer(maxParticleNum, Marshal.SizeOf(typeof(Uint2)));
            gridIndicesBuffer = new ComputeBuffer((int)(GridDim.x * GridDim.y), Marshal.SizeOf(typeof(Uint2)));
        }

        void PositionBasedFluid() {

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
            PBF_CS.SetVector("_GridDim", GridDim);
            PBF_CS.SetFloat("_GridH", gridH);


            int kernel = 0;

            // -----------------------------------------------------------------
            // Build Grid : 粒子の位置からグリッドハッシュとパーティクルIDを結びつける
            // -----------------------------------------------------------------
            kernel = PBF_CS.FindKernel("BuildGridCS");
            PBF_CS.SetBuffer(kernel, "_ParticlesBufferRead", particleBufferRead);
            PBF_CS.SetBuffer(kernel, "_GridBufferWrite", gridBuffer);
            PBF_CS.Dispatch(kernel, threadGroupSize, 1, 1);

            // -----------------------------------------------------------------
            // Sort Grid : グリッドインデックス順に粒子インデックスをソートする
            // -----------------------------------------------------------------
            GPUSort(gridBuffer, gridPingPongBuffer);    // TODO ソートされてない??

            // -----------------------------------------------------------------
            // Build Grid Indices : グリッドの開始終了インデックスを格納
            // -----------------------------------------------------------------
            // 初期化
            kernel = PBF_CS.FindKernel("ClearGridIndicesCS");
            PBF_CS.SetBuffer(kernel, "_GridIndicesBufferWrite", gridIndicesBuffer);
            PBF_CS.Dispatch(kernel, (int)(GridDim.x * GridDim.y) / SIMULATION_BLOCK_SIZE, 1, 1);

            // 格納
            kernel = PBF_CS.FindKernel("BuildGridIndicesCS");
            PBF_CS.SetBuffer(kernel, "_GridBufferRead", gridBuffer);
            PBF_CS.SetBuffer(kernel, "_GridIndicesBufferWrite", gridIndicesBuffer);
            PBF_CS.Dispatch(kernel, threadGroupSize, 1, 1);
            

            //　-----------------------------------------------------------------
            // Rearrange : ソートしたグリッド関連付け配列からパーティクルIDだけを取り出す
            //　-----------------------------------------------------------------

            kernel = PBF_CS.FindKernel("RearrangeParticlesCS");
            PBF_CS.SetBuffer(kernel, "_GridBufferRead", gridBuffer);
            PBF_CS.SetBuffer(kernel, "_ParticlesBufferRead", particleBufferRead);
            PBF_CS.SetBuffer(kernel, "_ParticlesBufferWrite", sortedParticleBuffer);
            PBF_CS.Dispatch(kernel, threadGroupSize, 1, 1);
            

            // Copy Buffer Every Frame
            kernel = PBF_CS.FindKernel("CopyBuffer");
            PBF_CS.SetBuffer(kernel, "_ParticlesBufferRead", sortedParticleBuffer);
            PBF_CS.SetBuffer(kernel, "_ParticlesBufferWrite", particleBufferWrite);
            PBF_CS.Dispatch(kernel, threadGroupSize, 1, 1);

            // Process1. Update
            kernel = PBF_CS.FindKernel("Update");
            PBF_CS.SetBuffer(kernel, "_ParticlesBufferRead", sortedParticleBuffer);
            PBF_CS.SetBuffer(kernel, "_ParticlesBufferWrite", particleBufferWrite);
            //PBF_CS.SetBuffer(kernel, "_GridIndicesBufferRead", gridIndicesBuffer);
            PBF_CS.Dispatch(kernel, threadGroupSize, 1, 1);


            // Process2 Lambda Iteration
            for (int i = 0; i < PRESSURE_ITERATIONS; i++) {

                // CalcLambda
                kernel = PBF_CS.FindKernel("CalcLambda");
                PBF_CS.SetBuffer(kernel, "_ParticlesBufferRead", sortedParticleBuffer);
                PBF_CS.SetBuffer(kernel, "_ParticlesBufferWrite", particleBufferWrite);
                PBF_CS.SetBuffer(kernel, "_GridIndicesBufferRead", gridIndicesBuffer);
                PBF_CS.Dispatch(kernel, threadGroupSize, 1, 1);

                // CalcDeltaP
                kernel = PBF_CS.FindKernel("CalcDeltaP");
                PBF_CS.SetBuffer(kernel, "_ParticlesBufferRead", sortedParticleBuffer);
                PBF_CS.SetBuffer(kernel, "_ParticlesBufferWrite", particleBufferWrite);
                PBF_CS.SetBuffer(kernel, "_GridIndicesBufferRead", gridIndicesBuffer);
                PBF_CS.Dispatch(kernel, threadGroupSize, 1, 1);

                // ApplyDeltaP
                kernel = PBF_CS.FindKernel("ApplyDeltaP");
                PBF_CS.SetBuffer(kernel, "_ParticlesBufferRead", sortedParticleBuffer);
                PBF_CS.SetBuffer(kernel, "_ParticlesBufferWrite", particleBufferWrite);
                PBF_CS.SetBuffer(kernel, "_GridIndicesBufferRead", gridIndicesBuffer);
                PBF_CS.Dispatch(kernel, threadGroupSize, 1, 1);
            }

            // Update Velocity
            kernel = PBF_CS.FindKernel("UpdateVelocity");
            PBF_CS.SetBuffer(kernel, "_ParticlesBufferRead", sortedParticleBuffer);
            PBF_CS.SetBuffer(kernel, "_ParticlesBufferWrite", particleBufferWrite);
            PBF_CS.SetBuffer(kernel, "_GridIndicesBufferRead", gridIndicesBuffer);
            PBF_CS.Dispatch(kernel, threadGroupSize, 1, 1);

            SwapBuffer(ref particleBufferRead, ref particleBufferWrite);
        }

        void CreateWater() {
            
            for (int i = 0; i < maxParticleNum; i++) {
                particles.Add(new Particle(new Vector2(Random.value * range.x / 2, Random.value * Mathf.Min(range.y, 256)), 1));
            }
            

            particle_array = particles.ToArray();
            particleBufferRead.SetData(particle_array);
            particleBufferWrite.SetData(particle_array);
            threadGroupSize = maxParticleNum / SIMULATION_BLOCK_SIZE;
        }

        void GPUSort(ComputeBuffer inBuffer, ComputeBuffer tempBuffer) {
            ComputeShader sortCS = BitonicCS;

            int KERNEL_ID_BITONICSORT = sortCS.FindKernel("BitonicSort");
            int KERNEL_ID_TRANSPOSE = sortCS.FindKernel("MatrixTranspose");

            uint NUM_ELEMENTS = (uint)maxParticleNum;
            uint MATRIX_WIDTH = BITONIC_BLOCK_SIZE;
            uint MATRIX_HEIGHT = (uint)NUM_ELEMENTS / BITONIC_BLOCK_SIZE;

            for (uint level = 2; level <= BITONIC_BLOCK_SIZE; level <<= 1) {
                SetGPUSortConstants(sortCS, level, level, MATRIX_HEIGHT, MATRIX_WIDTH);

                // Sort the row data
                sortCS.SetBuffer(KERNEL_ID_BITONICSORT, "Data", inBuffer);
                sortCS.Dispatch(KERNEL_ID_BITONICSORT, (int)(NUM_ELEMENTS / BITONIC_BLOCK_SIZE), 1, 1);
            }

            // Then sort the rows and columns for the levels > than the block size
            // Transpose. Sort the Columns. Transpose. Sort the Rows.
            for (uint level = (BITONIC_BLOCK_SIZE << 1); level <= NUM_ELEMENTS; level <<= 1) {
                // Transpose the data from buffer 1 into buffer 2
                SetGPUSortConstants(sortCS, level / BITONIC_BLOCK_SIZE, (level & ~NUM_ELEMENTS) / BITONIC_BLOCK_SIZE, MATRIX_WIDTH, MATRIX_HEIGHT);
                sortCS.SetBuffer(KERNEL_ID_TRANSPOSE, "Input", inBuffer);
                sortCS.SetBuffer(KERNEL_ID_TRANSPOSE, "Data", tempBuffer);
                sortCS.Dispatch(KERNEL_ID_TRANSPOSE, (int)(MATRIX_WIDTH / TRANSPOSE_BLOCK_SIZE), (int)(MATRIX_HEIGHT / TRANSPOSE_BLOCK_SIZE), 1);

                // Sort the transposed column data
                sortCS.SetBuffer(KERNEL_ID_BITONICSORT, "Data", tempBuffer);
                sortCS.Dispatch(KERNEL_ID_BITONICSORT, (int)(NUM_ELEMENTS / BITONIC_BLOCK_SIZE), 1, 1);

                // Transpose the data from buffer 2 back into buffer 1
                SetGPUSortConstants(sortCS, BITONIC_BLOCK_SIZE, level, MATRIX_HEIGHT, MATRIX_WIDTH);
                sortCS.SetBuffer(KERNEL_ID_TRANSPOSE, "Input", tempBuffer);
                sortCS.SetBuffer(KERNEL_ID_TRANSPOSE, "Data", inBuffer);
                sortCS.Dispatch(KERNEL_ID_TRANSPOSE, (int)(MATRIX_HEIGHT / TRANSPOSE_BLOCK_SIZE), (int)(MATRIX_WIDTH / TRANSPOSE_BLOCK_SIZE), 1);

                // Sort the row data
                sortCS.SetBuffer(KERNEL_ID_BITONICSORT, "Data", inBuffer);
                sortCS.Dispatch(KERNEL_ID_BITONICSORT, (int)(NUM_ELEMENTS / BITONIC_BLOCK_SIZE), 1, 1);
            }
        }

        void SetGPUSortConstants(ComputeShader cs, uint level, uint levelMask, uint width, uint height) {
            cs.SetInt("_Level", (int)level);
            cs.SetInt("_LevelMask", (int)levelMask);
            cs.SetInt("_Width", (int)width);
            cs.SetInt("_Height", (int)height);
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
        public Vector3 color;

        public Particle(Vector2 pos, float mass) {
            this.oldPos = pos;
            this.mass = mass;
            this.lambda = 0;
            this.pConstraint = 0;
            this.newPos = new Vector2(0f, 0f);
            this.velocity = new Vector2(0f, 0f);
            this.force = new Vector2(0f, 0f);
            this.deltaP = new Vector2(0f, 0f);
            this.color = new Vector3(1, 1, 1);
        }
    }
    
    struct Uint2 {
        public uint x;
        public uint y;
    }

}
