using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PBD : MonoBehaviour {

    private List<Particle> particles = new List<Particle>();
    private CellGrid grid;

    public Vector3 GRAVITY = new Vector3(0f, -9.8f, 0f);

    //Number of iterations to update pressure constraints (Macklin et al. used 4)
    public int PRESSURE_ITERATIONS = 6;
    public float dt = 0.005f;

    private static float H = 1.25f;
    private static float KPOLY = (float)(315f / (64f * Mathf.PI * Mathf.Pow(H, 9)));
    private static float SPIKY = (float)(45f / (Mathf.PI * Mathf.Pow(H, 6)));
    private static float VISC = (float)(15f / (2 * Mathf.PI * (H * H * H)));
    private static float REST_DENSITY = 1f;

    // Epsilon used in lambda calculation
    // See Macklin part 3
    private static float EPSILON_LAMBDA = 150f;
    private static float C = 0.01f;

    // K and deltaQMag used in sCorr Calculation
    // See Macklin part 4
    private static float EPSILON_VORTICITY = 10f;
    private static float K = 0.001f;
    private static float deltaQMag = 0.3f * H;
    private static float wQH = KPOLY * (H * H - deltaQMag * deltaQMag) * (H * H - deltaQMag * deltaQMag) * (H * H - deltaQMag * deltaQMag);

    // Used for bounds of the box
    public float rangex = 30f;
    public float rangey = 30f;
    public float rangez = 30f;

    public float time = 0;
    public bool dropped = false;

    public bool random_start = false;

    #region Rendering
    Mesh mesh;
    public Material displayMat;
    MaterialPropertyBlock _block;
    public MaterialPropertyBlock block {
        get {
            if (_block == null) {
                _block = new MaterialPropertyBlock();
            }
            return _block;
        }
    }
    #endregion Rendering


    void Start () {
        CreateWater();
        CreateGrid();
	}
	
	void Update () {

        PositionBasedDynamics();

        Render();
        
    }

    void OnDrawGizmos() {
        Gizmos.DrawWireCube(new Vector3(rangex, rangey, rangez)/2, new Vector3(rangex, rangey, rangez));
    }

    void Render() {
        CreateMesh();
        Graphics.DrawMesh(mesh, transform.localToWorldMatrix, displayMat, 0, null, 0, block);
    }

    void CreateMesh() {
        mesh = BuildMesh(GetPositions().Count);
    }

    Mesh BuildMesh(int num_particles) {
        Mesh particleMesh = new Mesh();

        var vertices = new Vector3[num_particles];
        var uvs = new Vector2[num_particles];
        var indices = new int[num_particles];

        for (int i = 0; i < num_particles; i++) {
            vertices[i] = particles[i].oldPos;
            uvs[i] = new Vector2(0, 0);
            indices[i] = i;
        }

        particleMesh.vertices = vertices;
        particleMesh.uv = uvs;

        particleMesh.SetIndices(indices, MeshTopology.Points, 0);

        return particleMesh;
    }

    void PositionBasedDynamics() {
        foreach (Particle p in particles) {
            ApplyGravity(p);
            p.setNewPos(p.getOldPos());

            // update velocity vi = vi + delta T * fext
            p.velocity += p.getForce() * dt;

            // predict position x* = xi + delta T * vi
            p.newPos += p.getVelocity() * dt;

            ImposeConstraints(p);
        }

        // get neighbors
        grid.updateCells(particles);

        foreach (Particle p in particles) {
            List<Particle> neighborParticles = new List<Particle>();
            List<Cell> neighborCells = p.getCell().getNeighbors();
            foreach (Cell c in neighborCells) {
                List<Particle> allParticles = c.getParticles();
                List<Particle> nearParticles = new List<Particle>();
                foreach (Particle n in allParticles) {
                    if (Vector3.Distance(p.getNewPos(), n.getNewPos()) <= H) {
                        nearParticles.Add(n);
                    }
                }
                neighborParticles.AddRange(nearParticles);  // ?? addAll
            }

            neighborParticles.Remove(p);
            p.setNeighbors(neighborParticles);
        }

        // while solver < iterations (they say that 2-4 is enough in the paper)
        for (int i = 0; i < PRESSURE_ITERATIONS; i++) {
            // Set lambda
            foreach (Particle p in particles) {
                List<Particle> neighbors = p.getNeighbors();
                p.setLambda(Lambda(p, neighbors));
            }
            // Calculate deltaP
            foreach (Particle p in particles) {
                Vector3 deltaP = new Vector3(0f, 0f, 0f);
                List<Particle> neighbors = p.getNeighbors();
                foreach (Particle n in neighbors) {
                    float lambdaSum = p.getLambda() + n.getLambda();
                    float sCorr = SCorr(p, n);
                    // float sCorr = 0;
                    deltaP += WSpiky(p.getNewPos(), n.getNewPos()) * (lambdaSum + sCorr);
                }

                p.setDeltaP(deltaP / REST_DENSITY);
            }

            // Update position x*i = x*i + delta Pi
            foreach (Particle p in particles) {
                p.newPos += p.getDeltaP();
            }
        }

        foreach (Particle p in particles) {
            ImposeConstraints(p);

            // set new velocity vi = (x*i - xi) /(deltaT) 
            p.setVelocity((p.getNewPos() - p.getOldPos()) / dt);

            // apply vorticity confinement
            p.velocity += VorticityForce(p) * dt;

            // apply XSPH viscosity
            p.velocity += XsphViscosity(p);

            // update position xi = x*i
            p.setOldPos(p.getNewPos());
        }
    }

    void CreateWater() {
        if (!random_start) {
            for (int i = 1; i < 30; i++) {
                for (int j = 10; j < 30; j++) {
                    for (int k = 1; k < 10; k++) {
                        particles.Add(new Particle(new Vector3(i, j, k), 1f));
                    }
                }
            }
        } else {
            for (int i = 0; i < 5000; i++) {
                particles.Add(new Particle(new Vector3(Random.value * (float)rangex, Random.value * (float)rangey, Random.value * (float)rangez), 1));
            }
        }
    }

    void CreateGrid() {
        // create cell grid
        grid = new CellGrid((int)rangex, (int)rangey, (int)rangez); // should be whatever the size of the box is
    }

    public List<Vector3> GetPositions() {
        List<Vector3> positions = new List<Vector3>();

        foreach(Particle part in particles) {
            Vector3 pos = part.getOldPos();
            positions.Add(new Vector3(pos.x, pos.y, pos.z));
        }

        return positions;
    }



    private void ApplyGravity(Particle p) {
        p.setForce(0f, 0f, 0f);
        p.force += GRAVITY;
    }

    // Poly6 Kernel
    private float WPoly6(Vector3 pi, Vector3 pj) {
        Vector3 r = pi - pj;
        float rLen = r.magnitude;   // .Length();
        if (rLen > H || rLen == 0) {
            return 0;
        }
        return (float)(KPOLY * Mathf.Pow((H * H - r.sqrMagnitude), 3));
    }

    // Spiky Kernel
    private Vector3 WSpiky(Vector3 pi, Vector3 pj) {
        Vector3 r = pi - pj;
        float rLen = r.magnitude;
        if (rLen > H || rLen == 0) {
            return new Vector3(0f, 0f, 0f);
        }

        float coeff = (H - rLen) * (H - rLen);
        coeff *= SPIKY;
        coeff /= rLen;
        return r * (-1f * coeff);
    }

    private Vector3 WViscosity(Vector3 pi, Vector3 pj) {
        Vector3 r = pi - pj;
        float rLen = r.magnitude;
        if (rLen > H || rLen == 0) return new Vector3(0f, 0f, 0f);

        float coeff = (-1 * (rLen * rLen * rLen)) / (2 * (H * H * H));
        coeff += (r.sqrMagnitude / (H * H));
        coeff += (H / (2 * rLen)) - 1;
        return r * coeff;
    }

    //Calculate the lambda value for pressure corrections
    private float Lambda(Particle p, List<Particle> neighbors) {
        float densityConstraint = CalcDensityConstraint(p, neighbors);
        Vector3 gradientI = new Vector3(0f, 0f, 0f);
        float sumGradients = 0;
        foreach (Particle n in neighbors) {
            // Calculate gradient with respect to j
            Vector3 gradientJ = WSpiky(p.getNewPos(), n.getNewPos()) / REST_DENSITY;

            // Add magnitude squared to sum
            sumGradients += gradientJ.sqrMagnitude;
            // Continue calculating particle i gradient
            gradientI += gradientJ; // ?? add
        }
        // Add the particle i gradient magnitude squared to sum
        sumGradients += gradientI.sqrMagnitude;
        return ((-1f) * densityConstraint) / (sumGradients + EPSILON_LAMBDA);
    }

    //Returns density constraint of a particle
    private float CalcDensityConstraint(Particle p, List<Particle> neighbors) {
        float sum = 0f;
        foreach (Particle n in neighbors) {
            sum += n.getMass() * WPoly6(p.getNewPos(), n.getNewPos());
        }

        return (sum / REST_DENSITY) - 1;
    }

    //Returns vorticity vector for a given particle
    private Vector3 Vorticity(Particle p) {
        Vector3 vorticity = new Vector3(0, 0, 0);
        Vector3 velocityDiff;
        Vector3 gradient;

        List<Particle> neighbors = p.getNeighbors();
        foreach (Particle n in neighbors) {
            velocityDiff = n.getVelocity() - p.getVelocity();
            gradient = WViscosity(p.getNewPos(), n.getNewPos());
            vorticity += Vector3.Cross(velocityDiff, gradient);  // add, CROSS?
        }

        return vorticity;
    }
    //Returns the eta vector that points in the direction of the corrective force
    private Vector3 Eta(Particle p, float vorticityMag) {
        List<Particle> neighbors = p.getNeighbors();
        Vector3 eta = new Vector3(0, 0, 0);
        foreach (Particle n in neighbors) {
            eta += WViscosity(p.getNewPos(), n.getNewPos()) * vorticityMag;
        }

        return eta;
    }

    //Calculates vorticity force for a particle
    private Vector3 VorticityForce(Particle p) {
        Vector3 vorticity = Vorticity(p);
        if (vorticity.magnitude == 0) {
            // No direction for eta
            return new Vector3(0f, 0f, 0f);
        }
        Vector3 eta = Eta(p, vorticity.magnitude);
        Vector3 n = eta.normalized;
        return Vector3.Cross(n, vorticity) * EPSILON_VORTICITY;
    }

    // Make sure that particle does not leave the cube grid
    private void ImposeConstraints(Particle p) {
        if (OutOfRange(p.getNewPos().x, 0, rangex)) {
            p.velocity.x = 0;
        }

        if (OutOfRange(p.getNewPos().y, 0, rangey)) {
            p.velocity.y = 0;
        }

        if (OutOfRange(p.getNewPos().z, 0, rangez)) {
            p.velocity.z = 0;
        }

        p.newPos.x = ClampedConstraint(p.getNewPos().x, rangex);
        p.newPos.y = ClampedConstraint(p.getNewPos().y, rangey);
        p.newPos.z = ClampedConstraint(p.getNewPos().z, rangez);
    }


    private float ClampedConstraint(float x, float max) {
        if (x < 0) {
            return 0;
        } else if (x >= max) {
            return max - 1e-3f;
        } else {
            return x;
        }
    }

    private float SCorr(Particle pi, Particle pj) {
        // Get Density from WPoly6 and divide by constant from paper
        float corr = WPoly6(pi.getNewPos(), pj.getNewPos()) / wQH;
        // take to power of 4
        corr *= corr * corr * corr;
        return -K * corr;
    }

    //Returns Vector3 to add to velocity
    //See Macklin Part 5, Equation 17
    private Vector3 XsphViscosity(Particle p) {
        Vector3 visc = Vector3.zero;
        List<Particle> neighbors = p.getNeighbors();
        foreach (Particle n in neighbors) {
            Vector3 velocityDiff = n.getVelocity() - p.getVelocity();
            velocityDiff *= WPoly6(p.getNewPos(), n.getNewPos());
        }
        return visc *= C;
    }

    //Tests if a particle is out of range of the box
    private static bool OutOfRange(float x, float min, float max) {
        return x < min || x >= max;
    }

    //private void DropParticles() {
    //    for (int i = 25; i < 34; i++) {
    //        for (int j = 15; j < 50; j++) {
    //            for (int k = 0; k < 5; k++) {
    //                particles.Add(new Particle(new Vector3(i, j, k), 1f));
    //            }
    //        }
    //    }
    //}
}

