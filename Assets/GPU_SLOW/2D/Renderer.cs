using UnityEngine;

namespace PBF_GPU_SLOW_2D {

    public class Renderer : MonoBehaviour {

        public PBF GPUScript;

        public Material ParticleRenderMat;

        void OnRenderObject(){
            DrawObject();
        }

        void DrawObject(){
            Material m = ParticleRenderMat;
            m.SetPass(0);
            m.SetBuffer("_Particles", GPUScript.GetBuffer());
            Graphics.DrawProcedural(MeshTopology.Points, GPUScript.GetMaxParticleNum());
        }

    }

}