using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class TubCreator : MonoBehaviour {

    public float width = 30;
    public float height = 5;

    public int num_layers = 1;

    public float radius = 0.5f;
    private float diameter;

    public GameObject sphere;

    public Vector3 pos;

	void Start () {
        diameter = radius * 2;
        pos = transform.position - new Vector3(width / 2, 0, 0);

        int num_height = (int)(height / diameter);
        int num_width = (int)(width / diameter);

        // 下の壁
        for(int l = 0; l < num_layers; l++) {
            for(int i = 0; i<num_width; i++) {
                // pos + vec2(i * diameter + radius, -radius)
                Instantiate(sphere, pos + new Vector3(i * diameter + radius, -diameter * l -radius, 0), Quaternion.identity, transform);
            }
        }

        // 左の壁
        for (int l = 0; l < num_layers; l++) {
            for (int i = 0; i < num_height; i++) {
                // pos + vec2(i * diameter + radius, -radius)
                Instantiate(sphere, pos + new Vector3(diameter * l + radius, i * diameter + radius, 0), Quaternion.identity, transform);
            }
        }

        // 右の壁
        for (int l = 0; l < num_layers; l++) {
            for (int i = 0; i < num_height; i++) {
                Vector3 pivot = pos + new Vector3((num_width-num_layers)*diameter, pos.y, 0);
                Instantiate(sphere, pivot + new Vector3(diameter * l + radius, i * diameter + radius, 0), Quaternion.identity, transform);
            }
        }
    }

    void OnDrawGizmos() {
        Gizmos.DrawLine(pos + new Vector3(-0.5f, 0, 0), pos + new Vector3(0.5f, 0, 0));
        Gizmos.DrawLine(pos + new Vector3(0, -0.5f, 0), pos + new Vector3(0, 0.5f, 0));

        Gizmos.color = Color.red;
        Gizmos.DrawLine(transform.position + new Vector3(-0.5f, 0, 0), transform.position + new Vector3(0.5f, 0, 0));
        Gizmos.DrawLine(transform.position + new Vector3(0, -0.5f, 0), transform.position + new Vector3(0, 0.5f, 0));
    }
	
	void Update () {
		
	}
}
