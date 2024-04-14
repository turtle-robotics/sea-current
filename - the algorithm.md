- the algorithm
    - loop through the obstacles
    - grab the end points
    - pick an end point and put the rest into a stack
    - using the picked point, derive "look-around" point from secant line
    # we now only care abt the one obstacle
    - put "blinders" on and trace the obstacle
    - if obstacle is closed pop the stack and trace another obstacle
    - otherwise we end up tracing the maze ggwp

- changes to the code
    - add closed/open data to obstacle
    - edit the "contains" function to add support for open obstacles
    - add endpoint data for open obstacles

Movment Algorithum
    - 3 states
        - State 1: "Lost Bug"   Inputs(NONE), Output(Move Forward)
            - If robot has no object data to move to it will move forward until an object is found by    the LIDAR
        - State 2: "Curious Bug"    Inputs(A Specific Target Point), Output()
            - Once the robot has object data it will default to this case. 
            - In this case the robot will move to a object target point to map the rest of the object. 
            - Once the robot has the target point it will create a straigt line called "m"
            - The robot will follow "m" to the target point, unless it encounters an object, to which it will then go around until it can resume its path on "m"
            - Once at this point the robot will then enter State 3
        - State 3: "Obsessed Bug"
            - The robot will first create a start point at its current position. 
            - It will then travel around the object until it has reached the the starting point it created.
            - Once the robot has reached the outer wall, this state will result in accomplishing the goal
