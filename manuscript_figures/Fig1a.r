library(rgdal)
library(plotrix)
library(rgeos)
library(raster)
library(data.table)
library(magrittr)

sites <- fread('../data/tanzania.samples.meta.csv') %>%
         .[, .(Lat = unique(latitude), Long = unique(longitude)), by = c('location', 'country')]
sites[, species := 'arabiensis']

pal <- c(arabiensis = rgb(0.05, 0.6, 0.05))

# Load the map of Africa
Afmap <- readOGR(dsn = '../sample_location_map/Africa_SHP')

sites[, ':='(label.offset.lat = -0.65,
             label.offset.long = c(-0.5, 0.4))
]

source('../shared_functions/R_plotting.r')

shadow.text <- function(x, y, label, shadow.size, shadow.col = 'white', text.col = 'black', ...){
	for (i in c(1,1,-1,-1)){
		for (j in c(1,-1,1,-1)){
			 text(x + shadow.size*i, y + shadow.size*j, label, col = shadow.col, ...)
		}
	}
	text(x, y, label, col = text.col, ...)
}

# Plot the map 
plot.collection.sites <- function(sites, palette){
	sites <- copy(sites)
	par(mar = c(1, 1, 0.2, 0.2), mgp = c(0.5,0,0), tcl = -0.2)
	countries <-  unique(sites$country) %>%
	              {setNames(sub("'", "`", .), .)}
	plot(Afmap[Afmap$NA2_DESCRI %in% countries,], border = 'black', col = 'white', lwd = 1.5, bg = rgb(0.482, 0.694, 0.741))
	# Make the other countries grey
	plot(Afmap[!(Afmap$NA2_DESCRI %in% countries),], border = 'grey60', col = 'grey80', add = T, bty = 'l')
	# Draw the main countries again so the borders are bold all the way
	plot(Afmap[Afmap$NA2_DESCRI %in% countries,], border = 'black', col = 'white', lwd = 1.5, add = T)
	# Label the countires
	for (country in names(countries)){
		centroid = gCentroid(Afmap[Afmap$NA2_DESCRI == countries[country], ])
		shadow.text(centroid$x, centroid$y+0.5, country, 0.02, text.col = 'grey20', cex = 1.3, font = 2)
	}
	# Add the sampling locations. 
	sites[, 
		{draw.circle(Long, Lat, radius = 0.2, col = lighten.col(palette[species], 0.5), nv = 1000)
	 	 shadow.text(Long+label.offset.long, Lat+label.offset.lat, text.col = palette[species], location, 0.02, cex = 1)},
		by = 'location'
	]
	# Add coordinates and white border box to hide axis line
	axis(1, cex.axis = 0.7)
	axis(2, cex.axis = 0.7, padj = -0.3)
	p = par('usr')
	rect(p[1], p[3], p[2], p[4], lwd = 2, xpd = NA, border = 'white')
}

pdf('Fig1a.pdf', width = 2.5, height = 3)
plot.collection.sites(sites, pal)
scalebar(d = 200, # distance in km
         xy = c(1.9,4.5),
         type = "bar", 
         divs = 2, 
         below = "km", 
         lonlat = TRUE,
         label = c(0,100,200), 
         adj=c(0, -1), 
         lwd = 2)
dev.off()

