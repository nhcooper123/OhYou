rm( list = ls() )

layout( matrix( c(1:3 ), ncol = 1))

# rob.par™ - my trademark plotting parameters
rob.par <- par(bty = "l",
		cex.lab = 1.1,
		tcl = 0.3,
		las = 1,
		pch = 1,
		bty = "l")
	
# ---- Plot Tippy Tree -----
alpha <- read.csv( "Data/ProfileDataTippyOU50.csv" )
alpha <- alpha[1:500,]

plot( alpha[,1], alpha[,2] - max( alpha[,2] ), xlim = c(0,5), ylim = c(-20, 0), type = "l", xlab = "", ylab = "", main = "Tippy" )

for( i in 2:11) {
	lines( alpha[,1], alpha[,i]- max( alpha[,i] ) )
	}


lines( c(0, 5), c(-1.92, -1.92), col = "red", lty = "22" ) 



# ---- Plot Rooty Tree ----

alpha <- read.csv( "Data/ProfileDataRootyOU50.csv" )
alpha <- alpha[1:500,]

plot( alpha[,1], alpha[,2] - max( alpha[,2] ), xlim = c(0,5), ylim = c(-20, 0), type = "l", xlab = "", ylab = "log Likelihood" , main = "Rooty" )

for( i in 2:11) {
	lines( alpha[,1], alpha[,i]- max( alpha[,i] ) )
	}
	
lines( c(0, 5), c(-1.92, -1.92), col = "red", lty = "22" ) 

# ---- Plot Yule Tree ----
alpha <- read.csv( "Data/ProfileDataYuleOU50.csv" )
alpha <- alpha[1:500,]

plot( alpha[,1], alpha[,2] - max( alpha[,2] ), xlim = c(0,5), ylim = c(-20, 0), type = "l", xlab = expression(alpha), ylab = "" , main = "Yule")

for( i in 2:11) {
	lines( alpha[,1], alpha[,i]- max( alpha[,i] ) )
	}

lines( c(0, 5), c(-1.92, -1.92), col = "red", lty = "22" ) 


