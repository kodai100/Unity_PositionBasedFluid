using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace PBF_GPU_FAST_2D {

    public class WaveController : MonoBehaviour {

        public PBF script;

        float initialRangeX;
        float time;
        public float speed;
        public float waveStrength = 0.1f;

        void Start() {
            initialRangeX = script.range.x;
            time = 0;
        }
        
        void Update() {
            time += Time.deltaTime;
            script.range.x = (initialRangeX - waveStrength) - waveStrength * Mathf.Sin(speed * (time-5f));
        }
    }

}