#!/bin/ruby

require 'csv'

raw_dt	= 0.0005
to_dt	= 0.01

if to_dt < raw_dt
	puts "dt error"
end

ave_num = (1.0/raw_dt) / (1.0/to_dt).to_int

start_time	= 0.0 #249.2
end_time	= start_time + 6.0

now_time	= 0.0
thrust_sum	= 0.0
step		= 1

raw_csv = CSV.read(ARGV[0], headers: false)

raw_csv.each do | data |
	raw_time	= data[0].to_f
	raw_thrust	= data[1].to_f

	if raw_time >= start_time and raw_time < end_time
		time	= raw_time - start_time
		thrust	= raw_thrust

		thrust_sum += thrust

		if step % ave_num == 0
			now_thrust = thrust_sum / ave_num
			puts "#{now_time.round(2)},#{now_thrust}"
			thrust_sum = 0.0
		end

		if step % ave_num == 1
			now_time = time
		end

		step += 1
	end
end
