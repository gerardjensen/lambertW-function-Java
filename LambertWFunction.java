import org.apache.commons.math3.complex.Complex;

public class LambertWFunction {
    public static void main(String[] args) {
        Complex z = new Complex(-0.02, 0);
        System.out.println(LambertW(z, 1));
    }

    //z * exp(z)
    public static Complex zexpz(Complex z) {
        return z.multiply(z.exp());
    }

    //The derivative of z * exp(z) = exp(z) + z * exp(z)
    public static Complex zexpz_d(Complex z) {
        return z.exp().add(zexpz(z));
    }

    //The second derivative of z * exp(z) = 2. * exp(z) + z * exp(z)
    public static Complex zexpz_dd(Complex z) {
        return z.exp().multiply(2).add(zexpz(z));
    }

    //Determine the initial point for the root finding
    public static Complex InitPoint(Complex z, int k) {
        final double pi = 3.14159265358979323846;
        final double e = 2.71828182845904523536;
        Complex I = new Complex(0, 1);
        Complex two_pi_k_I = I.multiply(2).multiply(pi).multiply(k);
        Complex ip = z.log().add(two_pi_k_I).subtract(z.log().add(two_pi_k_I).log()); // initial point coming from the general asymptotic approximation
        Complex p = z.multiply(e).add(1).multiply(2).sqrt(); // used when we are close to the branch cut around zero and when k=0,-1

        if (z.subtract(-Math.exp(-1)).abs() <= 1) { //we are close to the branch cut, the initial point must be chosen carefully
            if (k == 0)
                ip = p.add(-1).subtract(p.pow(2).multiply(1 / 3d)).add(p.pow(3).multiply(11d / 72));
            if (k == 1 && z.getImaginary() < 0)
                ip = p.negate().add(-1).subtract(p.pow(2).multiply(1 / 3d)).add(p.pow(3).multiply(11d / 72));
            if (k == -1 && z.getImaginary() > 0)
                ip = p.negate().add(-1).subtract(p.pow(2).multiply(1 / 3d)).add(p.pow(3).multiply(11d / 72));
        }

        if (k == 0 && z.subtract(0.5).abs() <= 0.5) // (1,1) Pade approximant for W(0,a)
            ip = z.multiply(7.061302897).add(0.1237166).multiply(0.35173371)
                    .divide(z.multiply(2).add(1).multiply(0.827184).add(2));
        if (k == -1 && z.subtract(0.5).abs() <= 0.5) // (1,1) Pade approximant for W(-1,a)
            ip = I.multiply(4.22096).add(2.2591588985)
                    .multiply(I.multiply(33.767687754).negate().add(-14.073271).multiply(z)
                            .subtract(I.multiply(19.071643).negate().add(12.7127).multiply(z.multiply(2).add(1))))
                    .divide(I.multiply(10.629721).negate().add(17.23103).multiply(z.multiply(2).add(1)).negate().add(2))
                    .negate();
        return ip;
    }

    public static Complex LambertW(Complex z, int k) {
        //For some particular z and k W(z,k) has simple value:
        if (z.getReal() == 0 && z.getImaginary() == 0)
            return k == 0 ? new Complex(0) : Complex.INF.negate();
        if (z.getReal() == -Math.exp(-1) && z.getImaginary() == 0 && (k == 0 || k == -1))
            return new Complex(-1);
        if (z.getReal() == Math.exp(1) && z.getImaginary() == 0 && k == 0)
            return new Complex(1);

        //Halley method begins
        Complex w = InitPoint(z, k); // intermediate values in the Halley method
        Complex wprev = InitPoint(z, k);
        final int maxiter = 30; // max number of iterations. This eliminates improbable infinite loops
        int iter = 0; // iteration counter
        double prec = 1e-30; // difference threshold between the last two iteration results (or the iter number of iterations is taken)
        do {
            wprev = w;
            Complex d = zexpz(w).subtract(z).multiply(zexpz_d(w)).multiply(2)
                    .divide(zexpz_d(w).pow(2).multiply(2).subtract(zexpz(w).subtract(z).multiply(zexpz_dd(w))));
            w = w.subtract(d);
            iter++;
        } while (w.subtract(wprev).abs() > prec && iter < maxiter);
        return w;

    }
}
