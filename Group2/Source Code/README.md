#Notes to compile the code:
(1) Go to "project setting" => "common properties", and then remove "quantlib.vcproj";
(2) Add QuantLib and Boost path to the include & library path settings;
(3) Fix the linking of IRSCrefitModel.lib & CppUniTestFramework.lib in "project setting" => "linker" => "general" => "additional lib dir";
(4) Re-build;
