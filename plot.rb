class GeneratePlots
  attr_accessor :prog_name, :plot_name
  def initialize(prog_name, plot_name)
    @plot_name = plot_name
    @prog_name = prog_name
    ###### COMPILE PROOGRAMM #########
    system("make #{@plot_name}")
    system("make")
  end

  def generate_test_data
    system('echo > in 100')
    system("./#{prog_name} < in > data.out")
    system("rm", 'in')
  end

  def plot_h
    system("./#{plot_name} h < data.out > data_h.out")
    gnuplot_file = File.open("plot.out", "w")
    gnuplot_file.write("set xlabel 't'\n")
    gnuplot_file.write("set ylabel 'h'\n")
    gnuplot_file.write("plot 'data_h.out' with lines\n")
    gnuplot_file.write("pause 3\n")
    gnuplot_file.close()
    system('gnuplot plot.out')
  end

  def plot_s
    system("./#{plot_name} s < data.out > data_s.out")
    gnuplot_file = File.open("plot.out", "w")
    gnuplot_file.write("set xlabel 't'\n")
    gnuplot_file.write("set ylabel 'sigma_e'\n")
    gnuplot_file.write("plot 'data_s.out' with lines\n")
    gnuplot_file.write("pause 3\n")
    gnuplot_file.close()
    system("gnuplot plot.out")
  end

  def plot_debug
    system("./#{plot_name} d > data_d.out")
    gnuplot_file = File.open("plot.out", "w")
    gnuplot_file.write("set xlabel 'x_0'\n")
    gnuplot_file.write("set ylabel 'h'\n")
    gnuplot_file.write("plot 'data_d.out' with lines\n")
    gnuplot_file.write("pause -1\n")
    gnuplot_file.close()
    system("gnuplot plot.out")
  end

  def clean
    system('rm plot.out data.out data_d.out data_h.out data_s.out')
  end
end

gp = GeneratePlots.new('membrane', 'plot')
gp.generate_test_data
# gp.plot_h
# gp.plot_s
gp.plot_debug
# gp.clean