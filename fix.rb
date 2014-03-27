data = readlines
data.each_with_index do |d, i|
  q = d.split(" ")
  puts "#{q[0].to_f*(10)} #{(q[1].to_f-0.36520)*0.5}"
end
