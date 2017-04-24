using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;

public class DebugBitonicSort : MonoBehaviour {

    Uint2[] a;
    public int width    = 512;
    public int height   = 512;
    Texture2D texture1;
    Texture2D texture2;

    public ComputeShader BitonicCS;
    ComputeBuffer bufferRead;
    ComputeBuffer pingpong;
    static uint BITONIC_BLOCK_SIZE = 512;
    static uint TRANSPOSE_BLOCK_SIZE = 16;

    // Use this for initialization
    void Start () {

        texture1 = new Texture2D(width, height, TextureFormat.ARGB32, false);
        texture1.filterMode = FilterMode.Point;
        texture2 = new Texture2D(width, height, TextureFormat.ARGB32, false);
        texture2.filterMode = FilterMode.Point;

        bufferRead = new ComputeBuffer(width * height, Marshal.SizeOf(typeof(Uint2)));
        pingpong = new ComputeBuffer(width * height, Marshal.SizeOf(typeof(Uint2)));

        a = new Uint2[width * height];
        for(int i = 0; i<width * height; i++) {
            uint tmp = (uint)Random.Range(0f, 255f);
            a[i].x = tmp;
            a[i].y = tmp;
        }

        bufferRead.SetData(a);

        ApplyTexture();
	}
	
	void Update () {
        GPUSort(bufferRead, pingpong);
        ApplyTexture();
    }

    void GPUSort(ComputeBuffer inBuffer, ComputeBuffer tempBuffer) {
        ComputeShader sortCS = BitonicCS;

        int KERNEL_ID_BITONICSORT = sortCS.FindKernel("BitonicSort");
        int KERNEL_ID_TRANSPOSE = sortCS.FindKernel("MatrixTranspose");

        uint NUM_ELEMENTS = (uint)(width * height);
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

    void OnGUI() {
        GUI.DrawTexture(new Rect(new Vector2(0, 0), new Vector2(texture1.width, texture1.height)), texture1);
        GUI.DrawTexture(new Rect(new Vector2(width, 0), new Vector2(texture2.width, texture2.height)), texture2);
    }

    void OnDestroy() {
        bufferRead.Release();
        pingpong.Release();
    }

    void ApplyTexture() {
        bufferRead.GetData(a);
        for (int i = 0; i < width * height; i++) {
            texture1.SetPixel(i % width, (height - 1) - i / width, new Color((float)a[i].x/256, (float)a[i].x / 256, (float)a[i].x / 256));
            texture2.SetPixel(i % width, (height - 1) - i / width, new Color((float)a[i].y/256, (float)a[i].y/256, (float)a[i].y/ 256));
        }
        texture1.Apply();
        texture2.Apply();
    }

    struct Uint2 {
        public uint x;
        public uint y;
    }
}
