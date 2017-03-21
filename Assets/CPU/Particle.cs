using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace CPU_PBD {

    public class Particle {
        public Vector3 oldPos;
        public Vector3 newPos;
        public Vector3 velocity;
        public Vector3 force;
        public Vector3 deltaP;
        public float mass;
        public float lambda;
        public float pConstraint;
        public List<Particle> neighbors;
        public Cell cell;

        public Particle(Vector3 pos, float mass) {
            this.oldPos = pos;
            this.mass = mass;
            this.newPos = new Vector3(0f, 0f, 0f);
            this.velocity = new Vector3(0f, 0f, 0f);
            this.force = new Vector3(0f, 0f, 0f);
            this.deltaP = new Vector3(0f, 0f, 0f);
            this.neighbors = new List<Particle>();
        }

        public Vector3 getOldPos() {
            return oldPos;
        }

        public Vector3 getNewPos() {
            return newPos;
        }

        public Vector3 getVelocity() {
            return velocity;
        }

        public void setVelocity(Vector3 velocity) {
            this.velocity = velocity;
        }

        public Vector3 getForce() {
            return force;
        }

        public void setForce(float x, float y, float z) {
            force.x = x;
            force.y = y;
            force.z = z;
        }

        public float getMass() {
            return mass;
        }

        public List<Particle> getNeighbors() {
            return neighbors;
        }

        public void setNeighbors(List<Particle> neighbors) {
            this.neighbors = neighbors;
        }

        public float getPConstraint() {
            return pConstraint;
        }

        public void setPConstraint(float f) {
            pConstraint = f;
        }

        public void setNewPos(Vector3 v) {
            newPos = v;
        }

        public void setOldPos(Vector3 v) {
            oldPos = v;
        }

        public Vector3 getDeltaP() {
            return deltaP;
        }

        public void setDeltaP(Vector3 deltaP) {
            this.deltaP = deltaP;
        }

        public float getLambda() {
            return lambda;
        }

        public void setLambda(float lambda) {
            this.lambda = lambda;
        }

        public Cell getCell() {
            return cell;
        }

        public void setCell(Cell cell) {
            this.cell = cell;
        }
    }

}