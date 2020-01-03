1. NonRTOptimalGCSRotate is the latest version. Rotate and straight walk are both validated in the file. Relative parameters are defined in the PhD Thesis, such as C1, C2. Chinese instructions are added in NonRTOptimalGCSRotate.cpp.

2. NonRTOptimalGCS and NonRTOptimal360 only validates the situation of straight walk. NonRTOptimalBCS is the method developed in the body coordinate system, and I'm not sure whether it is effective now.

3. Refer to the function "GetTimeOptimalGait" in three cpp .files to learn how to use the class.

4. Copy all files in RobotVIII_demo/Server to Robots/demo/demo_Plan, enter Robots/build to make, and finally enter Robots/build/bin to execute ./demo_Plan. When making the project, don't forget to add relative file names into the CMakeLists.txt.
