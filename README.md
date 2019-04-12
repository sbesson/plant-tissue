# README

This repository contains a MATLAB collection of objects and methods to analyse
the division plane of plant cells and plant tissues.

## Installation

- Clone or download a copy of this repository
- Add the folder and its subfolders to the MATLAB path

## Cell network

A cell network can be representated as a collection of vertices, edges and cells.

![cell network](celldescription.svg)

Vertices are points in the 2D space and can be coded as an `nx2` array.

Edges are nx4 arrays defined by several parameters:

- the first element is the index of the starting vertex
- the second element is the index of the ending vertex
- the third element is the oriented angle formed between the arc and the chord length joining the 2 vertices
- the fourth element is an optional parameter that contains the age of the vertex.

Cells are implemented as arrays of cells. Each cell contain a vector with the indexes to the edges.

## Example

The following example create a triangular cell:

    % Create a list of 2D vertices
    V=[0 0;
      1 0;
      0 1;];

    % Create a list of edges
    E = [1 2 0;
         2 3 pi/4;
         3 1 0;];

    % Create a cell array of cells
    C = {[1 2 3]};

    % Create the cellNetwork object
    N = cellNetwork(V,E,C);

    % Plot the cellNetwork object
    h = plot(N, 'o-r')
