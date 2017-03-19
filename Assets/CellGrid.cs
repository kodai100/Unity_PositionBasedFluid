using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CellGrid {
    private int width;
    private int height;
    private int depth;
    private Cell[,,] cells; //fixed size uniform grid of cells

	public CellGrid(int width, int height, int depth) {
        this.width = width;
        this.height = height;
        this.depth = depth;
        cells = new Cell[width,height,depth];
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                for (int k = 0; k < depth; k++) {
                    cells[i,j,k] = new Cell();
                }
            }
        }

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                for (int k = 0; k < depth; k++) {
                    for (int x = -1; x < 2; x++) {
                        for (int y = -1; y < 2; y++) {
                            for (int z = -1; z < 2; z++) {
                                if (i + x >= 0 && i + x < width && j + y >= 0 && j + y < height && k + z >= 0 && k + z < depth) {
                                    cells[i,j,k].addNeighbor(cells[i + x,j + y,k + z]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // 粒子位置のセルに、自身を追加する
    public void updateCells(List<Particle> particles) {
        clearCells();
        foreach (Particle p in particles) {
            Vector3 pos = p.getNewPos();
            //assuming indices are always valid because the box keeps the particles contained
            Cell cell = cells[(int)pos.x, (int)pos.y, (int)pos.z];
            cell.addParticle(p);
            p.setCell(cell);
        }
    }

    public void clearCells() {
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                for (int k = 0; k < depth; k++) {
                    cells[i, j, k].clearParticles();
                }
            }
        }
    }

    public int getWidth() {
        return width;
    }

    public int getHeight() {
        return height;
    }

    public int getDepth() {
        return depth;
    }
}