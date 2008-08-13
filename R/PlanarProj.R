`PlanarProj` <-
function (lon, lat) 
{
    lon0 <- mean(lon)
    lat0 <- mean(lat)
    x <- (lon - lon0) * pi/180
    y <- log(tan((pi/4 + lat/2) * pi/180))
    y <- list(x = x, y = y)
}
