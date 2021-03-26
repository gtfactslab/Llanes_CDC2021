# quadrotor_asif_mm
quadrotor_asif_mm is a ROS package that implements the main algorithm in [1] in the form of a case study for collision avoidance between two multirotors.

## Dependencies
```bash
mavros
mavros-extras
Gazebo
[Recommended] catkin-build-tools
```
## Installation
```bash
catkin build quadrotor_asif_mm
```

## Launch
```bash
roslaunch quadrotor_asif_mm quad_asif.launch
```

## References
[1] C. Llanes, M. Abate, S. Coogan, "Safety from Fast, In-the-Loop Reachability with Application to UAVs," submitted to 2021 60th IEEE Conference on Decision and Control (CDC), 2021.