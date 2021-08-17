public class complex {
	private double a;
	private double b;

	//4 different constructors, for ease of instantiating
	public complex() {
		a = 0.0;
		b = 0.0;
	}

	public complex (double x, double y) {
		a = x;
		b = y;
	}

	public complex (double val) {
		a = val;
		b = 0;
	}

	public complex (complex x) {
		this.a = x.a;
		this.b = x.b;
	}

	//equals method, checks if the components of two complex numbers are equal
	public boolean equals(complex y) {
		double epsilon = 0.0001;
		boolean equal = false;
		if ( ((this.a - y.a)/this.a) < epsilon && ((this.b - y.b)/this.b) < epsilon ) {
				equal = true;
			}
		return equal;
	}

	public boolean complexequals(Object o) {
		double epsilon = 0.0001;
		boolean equal = false;
		if (this.getClass() == o.getClass()) {
			complex w = (complex) o;
			if ( ((this.a - w.a)/this.a) < epsilon && ((this.b - w.b)/this.b) < epsilon ) {
				equal = true;
			}
		}
		return equal;
	}

	//returns a complex number in string format
	public static String toString(complex x) {
		return x.a + " + " + x.b + "i";
	}

	public String toString() {
		return a + " + " + b + "i";
	}

	//getter methods to return values of complex numbers
	public double getReal() {
		return a;
	}

	public double getIm() {
		return b;
	}

	//setter methods to set new values of complex numbers
	public double setReal(double val) {
		a = val;
		return a;
	}

	public double setIm(double val) {
		b = val;
		return b;
	}

	public static double arg(complex z) {
		double c = 0;
		if (z.b < 0) {
			c = (2 * Math.PI) + Math.atan2(z.b, z.a);
		} else if (z.b > 0) {
			c = Math.atan2(z.b, z.a);
		}
		return c;
	}

	//returns theta of (a,b)
	public double arg() {
		double c = 0;
		if (b < 0) {
			c = (2 * Math.PI) + Math.atan2(b, a);
		} else if (b > 0) {
			c = Math.atan2(b, a);
		}
		return c;
	}

	//returns magnitude of coordinates (a,b)
	public double mag() {
		double c = Math.sqrt( (a * a) + (b * b) );
		return c;
	}

	public static double mag(complex z) {
		return Math.sqrt( (z.a * z.a) + (z.b * z.b) );
	}

	public static complex conjugate(complex x) {
		complex z = new complex(x.a, -x.b);
		return z;
	}
	//returns conjugate of a complex number
	public complex conjugate() {
		complex z = new complex(a, -b);
		return z;
	}
	//adds two complex numbers x+y
	public static complex add(complex x, complex y) {
		double alpha = x.a + y.a; 
		double beta = x.b + y.b;
		complex z = new complex(alpha, beta);
		return z;
	}

	public complex add(complex y) {
		double alpha = this.a + y.a; 
		double beta = this.b + y.b;
		complex z = new complex(alpha, beta);
		return z;
	}
	//subtracts x-y complex numbers
	public static complex subtract(complex x, complex y) {
		double alpha = x.a - y.a; 
		double beta = x.b - y.b;
		complex z = new complex(alpha, beta);
		return z;
	}

	public complex subtract(complex y) {
		double alpha = a - y.a; 
		double beta = b - y.b;
		complex z = new complex(alpha, beta);
		return z;
	}

	//multiplies complex number * int n 
	public static complex multiply(complex x, double y) {
		complex z = new complex (y*x.a, y*x.b);
		return z;
	}

	public complex multiply(int y) {
		complex z = new complex (y*a, y*b);
		return z;
	}

	//multiples two complex numbers using FOIL
	public static complex multiply(complex x, complex y) {
		double alpha = (x.a * y.a) - (x.b * y.b);
		double beta = (x.a * y.b) + (x.b * y.a);
		complex z = new complex(alpha, beta);
		return z;
	}

	public complex multiply(complex y) {
		double alpha = (a * y.a) - (b * y.b);
		double beta = (a * y.b) + (b * y.a);
		complex z = new complex(alpha, beta);
		return z;
	}
	//uses methods above to calculate x/y, using the conjugate
	public static complex divide(complex x, complex y) {
		complex zed = new complex(0.0, 0.0);
		complex conj = new complex(y.a, -y.b);
		complex num = multiply(x, conj);
		complex den = multiply(y, conj);
		if (den.a == 0) {
			System.out.println("Error: division by 0");
			zed.setReal(999.999);
			zed.setIm(999.999);
		} else {
			double alpha = (num.a/den.a);
			double beta = (num.b/den.a);
			zed.setReal(alpha);
			zed.setIm(beta);
		}
		return zed;
	}

	public complex divide(complex y) {
		complex z = new complex(0.0, 0.0);
		complex conj = new complex(y.a, -y.b);
		complex num = multiply(this, conj);
		complex den = multiply(y, conj);
		if (den.a == 0) {
			System.out.println("Error: division by 0");
			z.setReal(999.999);
			z.setReal(999.999);
		} else {
			double alpha = (num.a/den.a);
			double beta = (num.b/den.a);
			z.setReal(alpha);
			z.setIm(beta);
		}
		return z;
	}

	//demoivre theorem to calculated x^n
	public static complex power (complex x, double n) {
		double theta = arg(x);
		double r = mag(x);
		double alpha = (Math.pow(r,n)*Math.cos(n*theta));
		double beta = (Math.pow(r,n)*Math.sin(n*theta));
		complex z = new complex(alpha, beta);
		return z;

	}

	public complex power (int n) {
		double theta = this.arg();
		double r = this.mag();
		double alpha = (Math.pow(r,n)*Math.cos(n*theta));
		double beta = (Math.pow(r,n)*Math.sin(n*theta));
		complex z = new complex(alpha, beta);
		return z;

	}
	//uses a fairly simple formula to calculate complex x^y
	public static complex complexpow (complex x, complex y) {
		double r = mag(x);
		double theta = arg(x);
		double logR = Math.log(r);
		double e = Math.E;
		double f = (y.a * logR)  - (y.b * theta);
		double ef = Math.pow(e, f);
		double g = (y.a * theta) + (y.b * logR); 
		double alpha = ef * Math.cos(g);
		double beta  = ef * Math.sin(g);
		complex z = new complex(alpha, beta);
		return z;
	}

	public complex complexpow (complex y) {
		double r = mag(this);
		double theta = arg(this);
		double logR = Math.log(r);
		double e = Math.E;
		double f = (y.a * logR)  - (y.b * theta);
		double ef = Math.pow(e, f);
		double g = (y.a * theta) + (y.b * logR); 
		double alpha = ef * Math.cos(g);
		double beta  = ef * Math.sin(g);
		complex z = new complex(alpha, beta);
		return z;
	}
	//simple formula that calculates natural log of a complex number
	public static complex complexln (complex x) {
		double r = mag(x);
		double theta = arg(x);
		double alpha = Math.log(r);
		double beta = theta;
		complex z = new complex(alpha, beta);
		return z;
	}

	public complex complexln () {
		double r = mag(this);
		double theta = arg(this);
		double alpha = Math.log(r);
		double beta = theta;
		complex z = new complex(alpha, beta);
		return z;
	}

	/*computes the roots of a quadratic equation with any coefficient being complex or real. It basically works by doing the positive and negative versions of the quadratic formula for each root you get 
	from square rooting the discriminant, then feeds those into two arrays. the arrays are then merged together into "roots' which is returned as the answer.*/
	public static complex[] quadratic (complex a, complex b, complex c) {
		complex det = subtract(power(b, 2), (multiply(multiply(a,4),c)));
		complex[] sqrdet = kroot(det, 2);
		complex[] proots = new complex[sqrdet.length];
		complex[] nroots = new complex[sqrdet.length];
		complex[] roots = new complex[sqrdet.length * 2];
		complex negB = new complex(-b.a, -b.b);
		complex root;

		for (int i = 0; i < proots.length; i++) {
			root = divide(add(negB, sqrdet[i]), multiply(a, 2));
			proots[i] = root;
		} 

		for (int i = 0; i < nroots.length; i++) {
			root = divide(subtract(negB, sqrdet[i]), multiply(a, 2));
			nroots[i] = root;
		}

		int count = 0;
		for (int i = 0; i < proots.length; i++) {
			roots[i] = proots[i];
			count++;
		} 

		for (int i = 0; i < nroots.length; i++) {
			roots[count++] = proots[i];
		} 

		return roots;
	}

	public static complex[] kroot (complex x, int n) {
		double newN = n;
		double theta = arg(x);
		double r = mag(x);
		double newR = Math.pow(r, 1/newN);
		double alpha;
		double beta;
		double nTheta;
		complex[] roots = new complex[n];

		for (int k = 0; k < n; k++) {
			nTheta = (theta + (2 * Math.PI * k)) / newN;
			alpha = newR * Math.cos(nTheta);
			beta = newR * Math.sin(nTheta);
			roots[k] = new complex(alpha, beta);
		}

		return roots;	
	}

	//compute nth root of a complex number4. uses polar form and iterates theta + 2pik to get all the roots
	public complex[] kroot (int n) {
		double newN = n;
		double theta = arg(this);
		double r = mag(this);
		double newR = Math.pow(r, 1/newN);
		double alpha;
		double beta;
		double nTheta;
		complex[] roots = new complex[n];

		for (int k = 0; k < n; k++) {
			nTheta = (theta + (2 * Math.PI * k)) / newN;
			alpha = newR * Math.cos(nTheta);
			beta = newR * Math.sin(nTheta);
			roots[k] = new complex(alpha, beta);
		}

		return roots;	
	}

	//calculate cosine of a complex number
	public static complex complexcos (complex x) {
		double e = Math.E;

		//splitting formula into 4 "parts"- makes it easier to read and debug
		double part1 = (Math.pow(e,-x.b) * Math.cos(x.a));
		double part2 = (Math.pow(e,x.b) * Math.cos(x.a));
		double alpha = 0.5 * (part1 + part2);

		double part3 = (Math.pow(e,-x.b) * Math.sin(x.a));
		double part4 = (Math.pow(e,x.b) * -Math.sin(x.a));
		double beta  = 0.5 * (part3 + part4);

		complex z = new complex(alpha, beta);
		return z;
	}

	public complex complexcos () {
		double e = Math.E;

		double part1 = (Math.pow(e,-b) * Math.cos(a));
		double part2 = (Math.pow(e,b) * Math.cos(a));
		double alpha = 0.5 * (part1 + part2);

		double part3 = (Math.pow(e,-b) * Math.sin(a));
		double part4 = (Math.pow(e,b) * -Math.sin(a));
		double beta  = 0.5 * (part3 + part4);

		complex z = new complex(alpha, beta);
		return z;
	}

	public static complex complexsin (complex x) {
		double e = Math.E;

		double part1 = (Math.pow(e, -x.b) * Math.sin(x.a));
		double part2 = (Math.pow(e, x.b) * -Math.sin(x.a));
		double alpha = 0.5 * (part1 - part2);

		double part3 = (Math.pow(e, -x.b) * Math.cos(x.a));
		double part4 = (Math.pow(e, x.b) * Math.cos(x.a));;
		double beta = 0.5 * (part3 - part4);

		complex z = new complex(alpha, beta);
		return z;
	}

	//calculate sine of a complex number
	public complex complexsin () {
		double e = Math.E; //ease of coding

		double part1 = (Math.pow(e, -b) * Math.sin(a));
		double part2 = (Math.pow(e, b) * -Math.sin(a));
		//^^ splitting up the formula into 4 "parts"- makes it easier to read and debug
		double alpha = 0.5 * (part1 - part2);

		double part3 = (Math.pow(e, -b) * Math.cos(a));
		double part4 = (Math.pow(e, b) * Math.cos(a));;
		double beta = 0.5 * (part3 - part4);

		complex z = new complex(alpha, beta);
		return z;
	}

	public static complex[] complexarccos (complex x) {
		complex xsqr = power(x,2);
		complex sumX;
		complex sumX2;
		complex part1 = new complex(xsqr.a-1, xsqr.b); //z^2-1
		complex[] xroots = kroot(part1, 2); //sqrt(z^2-1)
		complex sumLN;
		complex sumLN2;
		complex[] answers = new complex[2];
		//now, we must add x+xroots = sumX -> full formula is -iln(z+sqrt(z^2-1))
		//then-iln(sumX), and feed them into our solutions "sols" array
		sumX = add(x, xroots[0]);
		sumLN = complexln(sumX);
		answers[0] = new complex(sumLN.b, -sumLN.a); //this is multiplying the * -i
		sumX2 = add(x, xroots[1]);
		sumLN2 = complexln(sumX2);
		answers[1] = new complex(sumLN2.b, -sumLN2.a);

		return answers;
	}

	public complex[] complexarccos () {
		complex xsqr = power(this,2);
		complex sumX;
		complex sumX2;
		complex part1 = new complex(xsqr.a-1, xsqr.b); //z^2-1
		complex[] xroots = kroot(part1, 2); //sqrt(z^2-1)
		complex sumLN;
		complex sumLN2;
		complex[] answers = new complex[2];
		//now, we must add x+xroots = sumX -> full formula is -iln(z+sqrt(z^2-1))
		//then-iln(sumX), and feed them into our solutions "sols" array
		sumX = add(this, xroots[0]);
		sumLN = complexln(sumX);
		answers[0] = new complex(sumLN.b, -sumLN.a); //this is multiplying the * -i
		sumX2 = add(this, xroots[1]);
		sumLN2 = complexln(sumX2);
		answers[1] = new complex(sumLN2.b, -sumLN2.a);

		return answers;
	}

	public static complex[] complexarcsin (complex x) {
		complex xsqr = power(x,2);
		complex sumX;
		complex sumX2;
		complex part1 = new complex(1-xsqr.a, -xsqr.b); //1-z^2
		complex[] xroots = kroot(part1, 2); //sqrt(1-z^2)
		complex sumLN;
		complex sumLN2;
		complex[] answers = new complex[2];
		//now, we must add ix+xroots = sumX -> full formula is -iln(iz+sqrt(1-z^2))
		//then-iln(sumX), and feed them into our solutions "sols" array
		complex ix = new complex(-x.b, x.a);
		sumX = add(ix, xroots[0]);
		sumLN = complexln(sumX);
		answers[0] = new complex(sumLN.b, -sumLN.a); //this is multiplying the * -i
		sumX2 = add(ix, xroots[1]);
		sumLN2 = complexln(sumX2);
		answers[1] = new complex(sumLN2.b, -sumLN2.a);

		return answers;
	}

	public complex[] complexarcsin () {
		complex xsqr = power(this,2);
		complex sumX;
		complex sumX2;
		complex part1 = new complex(1-xsqr.a, -xsqr.b); //1-z^2
		complex[] xroots = kroot(part1, 2); //sqrt(1-z^2)
		complex sumLN;
		complex sumLN2;
		complex[] answers = new complex[2];
		//now, we must add ix+xroots = sumX -> full formula is -iln(iz+sqrt(1-z^2))
		//then-iln(sumX), and feed them into our solutions "sols" array
		complex ix = new complex(-b, a);
		sumX = add(ix, xroots[0]);
		sumLN = complexln(sumX);
		answers[0] = new complex(sumLN.b, -sumLN.a); //this is multiplying the * -i
		sumX2 = add(ix, xroots[1]);
		sumLN2 = complexln(sumX2);
		answers[1] = new complex(sumLN2.b, -sumLN2.a);

		return answers;
	}
}

//arcsin, quadratic, nonstatic methods, documentation