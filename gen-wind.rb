#!/bin/ruby

require "csv"

data = CSV.read(ARGV[0], col_sep: " ", headers: false)

press_max = data[0][0].to_f
press_min = data[data.size-1][0].to_f

height = 0
press = 1013.25

while press > press_max
  r0=6356766
  geo_h = (r0*height*0.001) / (r0 + height*0.001)
  temperature = 15 - 6.5*geo_h + 273.15
  press = 1013.25 * ((288.15 / temperature)**(-5.256))
  #print "#{height} #{press}\n"
  height = height + 1
end

STDERR.puts "min height: #{height}"

i = 0
press_base = press_max
press_width = data[1][0].to_f - press_base
x_base = data[0][1].to_f
y_base = data[0][2].to_f
x_diff = data[1][1].to_f - x_base
y_diff = data[1][2].to_f - y_base
x = 0.0
y = 0.0

while true
  r0=6356766
  geo_h = (r0*height*0.001) / (r0 + height*0.001)
  temperature = 15 - 6.5*geo_h + 273.15
  press = 1013.25 * ((288.15 / temperature)**(-5.256))

  if press < press_min
    break
  elsif press <= data[i+1][0].to_f
    i = i + 1
    press_base = data[i][0].to_f
    x_base = data[i][1].to_f
    y_base = data[i][2].to_f
    x_diff = data[i+1][1].to_f - x_base
    y_diff = data[i+1][2].to_f - y_base
    press_width = data[i+1][0].to_f - press_base
    #STDERR.puts "#{x_base} #{x_diff}"
  end
  pdiff = press - press_base
  x = x_base + x_diff * (pdiff / press_width)
  y = y_base + y_diff * (pdiff / press_width)
  puts "#{x},#{y}"
  height = height + 1
end
