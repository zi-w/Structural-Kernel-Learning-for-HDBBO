#!/usr/bin/env python
# Copyright (c) 2017 Zi Wang
from push_world import *
import sys

def easy_push(rx, ry, init_angle, simu_steps, gx, gy):
    # the object initial position is at (0,0)
    # rx, ry is the initial position of the robot.
    # the robot has velocity (-rx, -ry), which is then normalized to 10
    # init angle is the angle of the robot
    # simu_steps is the number of steps simulated in box2d
    # (gx, gy) is the goal location of the object

    # robot initial position cannot be at (0,0), which is the object's location
    assert(rx != 0 or ry !=0)

    # Set the parameter to False if don't need gui
    world = b2WorldInterface(True)

    # object properties and robot properties
    oshape, osize, ofriction, odensity, bfriction, robot_shape, robot_size  = 'circle', 1, 0.01, 0.05, 0.01, 'rectangle', (0.3,1) 
    # you can change the initial location of the object
    thing_init_location = (0,0)

    thing,base = make_thing(500, 500, world, oshape, osize, ofriction, odensity, bfriction, thing_init_location)

    # robot velocity
    xvel = -rx;
    yvel = -ry;
    regu = np.linalg.norm([xvel,yvel])
    xvel = xvel / regu * 10;
    yvel = yvel / regu * 10;

    robot = end_effector(world, (rx,ry), base, init_angle, robot_shape, robot_size)

    # run simulation. It returns the object location
    ret = simu_push2(world, thing, robot, base, xvel, yvel, simu_steps)

    # get distance to the goal
    ret = np.linalg.norm(np.array([gx, gy]) - ret)

if __name__ == '__main__':
    rx = float(sys.argv[1])
    ry = float(sys.argv[2])
    init_angle = float(sys.argv[3])
    simu_steps = int(float(sys.argv[4]) * 10)
    gx = float(sys.argv[5])
    gy = float(sys.argv[6])
    # to make the push even easier, can set init_angle = np.arctan(ry/rx) if rx \neq 0
    ret = easy_push(rx, ry, init_angle, simu_steps, gx, gy)
    sys.stdout.write(str(ret))