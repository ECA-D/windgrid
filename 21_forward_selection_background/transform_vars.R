dmax = 150e3
station$d2c[(station$d2c>dmax)] = dmax
grid$D2C[(grid$D2C>dmax)] = dmax
station$d2c = sqrt(station$d2c)
grid$D2C = sqrt(grid$D2C)

station$slope = log(station$slope+0.001)
grid$SLOPE = log(grid$SLOPE+0.001)

station$avg_roughness = log(10/(station$avg_roughness+0.001))
grid$ROUGH = log(10/(grid$ROUGH+0.001))

station$era5_avg_windspeed = log(station$era5_avg_windspeed)
grid$ERA5SPEED = log(grid$ERA5SPEED)

station$avg_windspeed = log(station$avg_windspeed)

