data = readlines
data.each_with_index do |d, i|
  q = d.split(" ")
  puts "#{q[0].to_f*10*(i+1)} #{q[1].to_f-0.38888881}"
end
