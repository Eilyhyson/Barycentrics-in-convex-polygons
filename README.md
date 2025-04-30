

# Barycentric Coordinates: A Comprehensive Implementation

This repository is the implementation of my undergraduate thesis project on barycentric coordinates. The project originally stems from Dmitry Anisimov's work on Danston, but I've extended it by integrating several additional coordinate algorithms. These include centroid coordinates, maximum entropy coordinates, harmonic coordinates, local coordinates, PIC coordinates, and Napoleon triangle-derived MVC.

Additionally, the project calculates triangle subdivision faces and outputs them in a `.faces` file format.

## Features

- **Coordinate Algorithms**:
  - Mean value coordinates
  - Maximum Entropy Coordinates
  - Harmonic Coordinates
  - Local Coordinates
  - PIC Coordinates
  - MVC derived from Napoleon’s Triangle
- **Triangle Subdivision**:
  - Computes faces of triangle subdivisions and saves them in `.faces` files.

## Installation

### On Linux (Ubuntu)

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/barycentric-coordinates.git
   cd barycentric-coordinates
   ```

2. Ensure you have the necessary dependencies installed:
   ```bash
   sudo apt update
   sudo apt install build-essential cmake libboost-all-dev
   ```

3. Compile the project:
   ```bash
   mkdir build
   cd build
   cmake ..
   make
   ```

4. Run the project:
   ```bash
   ./barycentric_coordinates
   ```

### On macOS

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/barycentric-coordinates.git
   cd barycentric-coordinates
   ```

2. Install dependencies using Homebrew:
   ```bash
   brew install boost cmake
   ```

3. Compile the project:
   ```bash
   mkdir build
   cd build
   cmake ..
   make
   ```

4. Run the project:
   ```bash
   ./barycentric_coordinates
   ```

## Usage

After building the project, you can run the program to calculate the barycentric coordinates for a given set of points within a triangle. The output will include the barycentric coordinates for each point, as well as the faces of triangle subdivisions.

## Contributing

Feel free to fork this repository, submit pull requests, or raise issues if you have any questions or suggestions.

## License

This project is licensed under the MIT License.

