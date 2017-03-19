using System.Collections.Generic;

public class Cell {
    public List<Particle> particles;
    public List<Cell> neighbors;

    public Cell() {
        particles = new List<Particle>();
        neighbors = new List<Cell>();
    }

    public void addNeighbor(Cell cell) {
        neighbors.Add(cell);
    }

    public List<Cell> getNeighbors() {
        return neighbors;
    }

    public void setNeighbors(List<Cell> neighbors) {
        this.neighbors = neighbors;
    }

    public void addParticle(Particle particle) {
        particles.Add(particle);
    }

    public void clearParticles() {
        particles.Clear();
    }

    public List<Particle> getParticles() {
        return particles;
    }

    public void setParticles(List<Particle> particles) {
        this.particles = particles;
    }
}