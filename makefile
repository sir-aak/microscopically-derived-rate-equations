CC		= c++											# Compiler
CFLAGS	= -std=c++17 -Wall -Ofast -llapack -lopenblas	# Optionen
INCLUDE = -I /home/krimlowski/Andrej/armadillo-9.850.1/include/

all:	mdre

mdre:	mdre_main.o mdre_parameters.o mdre_command_line_options.o mdre_node_equations.o mdre_network.o mdre_solver.o mdre_time_series.o mdre_auxiliary.o mdre_analyzer.o mdre_TISI_finder.o mdre_features.o mdre_function_switches.o mdre_parameters.hpp mdre_command_line_options.hpp mdre_node_equations.hpp mdre_network.hpp mdre_solver.hpp mdre_time_series.hpp mdre_auxiliary.hpp mdre_analyzer.hpp mdre_TISI_finder.hpp mdre_features.hpp mdre_function_switches.hpp
	$(CC) $(CFLAGS) $(INCLUDE) -o mdre mdre_main.o mdre_parameters.o mdre_command_line_options.o mdre_node_equations.o mdre_network.o mdre_solver.o mdre_time_series.o mdre_auxiliary.o mdre_analyzer.o mdre_TISI_finder.o mdre_features.o mdre_function_switches.o

clean:
	rm -f mdre
	rm -f *.o

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c $<
