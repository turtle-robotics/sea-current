# Sea-Current
A path planning library for the the Texas A&M University Robotics Team and Leadership Experience Maze Solving Robot Subteam

## Algorithm
- A path consisting of joined line-segments is generated through the “Fast Matching Trees” algorithm
- For each pair of line segments, a bezier curve is generated to interpolate between the endpoints of the pair
- The bezier curve control points are adjusted to avoid any obstacles
- A velocity profile is generated to accompany the path

## Dependencies
- Eigen
- toppra (packaged in tree)
- nlohmann/json (packaged in tree)
- matplotlibcpp (for testing only; packaged in tree)
