using PackageCompiler
create_sysimage([:Pluto, :PlutoUI, :PlotlyLight, :LaTeXStrings, :Plots];
                # sysimage_path="sys_plots.so",
                sysimage_path="/home/jovyan/sysimage.so",
                #precompile_execution_file = "warmup.jl",
                #replace_default = true,
                cpu_target = PackageCompiler.default_app_cpu_target())
