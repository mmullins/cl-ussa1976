;; -*- fill-column:102; comment-column:24; -*-

;; US Standard Atmosphere, 1976

(in-package :cl-ussa1976)

(declaim (optimize (debug 3)))
;;(declaim (optimize (speed 3) (safety 0)))


;;;; NOTES

;; 'Standard' here refers to 'US Standard Atmosphere, 1976', NOAA-S/T 76-1562.

;; Values related to moles are per kmole, not (conventional) mole.

;; Constants +k+, +fNA+ and +Rstar+ (Boltzmann's constant, Avogadro's number and the universal gas
;; constant) are not quite equal to current CODATA values, but are those given in the standard.

;; Functions with '-lo' are applicable for Z <= 86 km (H <= 84.8520 km),
;; Functions with '-hi' are applicable for Z > 86 km.


;;;; UTILITIES

(defun interpolate-find-span (v x)
  "Return span, 0..(- (length v) 2) of vector V containing value X.
X is in span i if X >= (aref V i) and  X < (aref V (+ 1 i)).
If X < (aref v 0), 0 is returned.
If X > (aref v (- (length v) 2)), (- (length v) 2) is returned."
  (let ((n-2 (- (length v) 2)))
    (cond ((>= x (aref v n-2)) n-2)
          ((< x (aref v 1)) 0)
          (t (loop
                with lo fixnum = 0 ;?? lo = 1 ??
                with hi fixnum = n-2
                for mid fixnum = (floor (+ lo hi) 2) do
                  (cond ((< x (aref v mid)) (setf hi mid))
                        ((>= x (aref v (+ mid 1))) (setf lo mid))
                        (t (return mid))))))))

(defun interpolate (v-x v-y x)
  "Interpolate the data in the vectors V-X and V-Y to return the y value corresponding
to x value X."
  (let* ((span (interpolate-find-span v-x x))
         (x1 (aref v-x span))
         (x2 (aref v-x (+ span 1)))
         (a2 (/ (- x x1) (- x2 x1)))
         (a1 (- 1.0d0 a2)))
    (+ (* a1 (aref v-y span)) (* a2 (aref v-y (+ span 1))))))

(defun interpolate-uniform (x-start dx v-y x)
  "Interpolate the data in vector V-Y which corresponds to x values starting at X-START and uniformly
spaced by DX, to return the y value corresponding to x value X.  If X is beyond the ends of the x
range, return an extrapolation based on the appropriate end span."
  (let* ((max-span (- (length v-y) 2))
         (span (max 0 (min max-span (truncate (/ (- x x-start) dx)))))
         (x1 (+ x-start (* span dx)))
         (a2 (/ (- x x1) dx))
         (a1 (- 1.0d0 a2)))
    (+ (* a1 (aref v-y span)) (* a2 (aref v-y (+ span 1))))))

(defun interpolate-uniform-weights (x-start dx v-y x)
  "Return (values span a1 a2) for interpolation in vector V-Y corresponding to the location of X in
the x range starting at X-START with spans of length DX.  Span is 0..(- (length V-Y) 2).  a1 is the
weight to be applied to the value at index value span, a2 is the weight to be applied to the value at
span+1."
  (let* ((max-span (- (length v-y) 2))
         (span (max 0 (min max-span (truncate (/ (- x x-start) dx)))))
         (x1 (+ x-start (* span dx)))
         (a2 (/ (- x x1) dx))
         (a1 (- 1.0d0 a2)))
    (values span a1 a2)))

(defun trapezoid-uniform (dx v-y &key (start 0) (end (1- (length v-y))) (integral-start 0.0d0))
  "Return a vector which is the integral of V-Y vs uniformly spaced x values, spaced at DX, using the
trapezoidal rule.  If START and/or END are supplied, use them as limits in the input array V-Y over
which the integral is done.  Normally END would be greater than START, and the integration is from low
to high indexes.  However, if END is less than START, the integration starts at the higher index,
START, and proceeds to the lower index, END.  In either case, INTEGRAL-START, if supplied, is the
integration value at START.  Also, in either case, the sign of DX should correspond to increasing
indexes of V-Y.  If START and END define a range that is only a subset of the total range of V-Y, the
vector returned has just the values of the integral between and including START and END.  In that case
it will be smaller than V-Y."
  (let* ((v-integral (make-array (1+ (abs (- end start)))
                                 :element-type (array-element-type v-y))))
    (if (< end start)
        (loop for i from start downto end
           for integral = integral-start then (- integral
                                                 (* 0.5d0 dx (+ (aref v-y (1+ i)) (aref v-y i))))
           do (setf (aref v-integral i) integral))
        (loop for i from start to end
           for integral = integral-start then (+ integral
                                                 (* 0.5d0 dx (+ (aref v-y (1- i)) (aref v-y i))))
           do (setf (aref v-integral i) integral)))
    v-integral))


;;;; CONSTANTS

;; using defparameter because defconstant gives warning on every recompile which is a pain during
;; development; would probably be more efficient to change to defconstant

;; category I constants
(defparameter +k+     1.380622d-23     "Boltzmann's constant [Nm/K]")
(defparameter +M_N2+  28.0134d0        "molecular weight of N_2")
(defparameter +M_O2+  31.9988d0        "molecular weight of O_2")
(defparameter +M_Ar+  39.948d0         "molecular weight of Ar")
(defparameter +M_CO2+ 44.00995d0       "molecular weight of CO_2")
(defparameter +M_Ne+  20.183d0         "molecular weight of Ne")
(defparameter +M_He+  4.0026d0         "molecular weight of He")
(defparameter +M_Kr+  83.80d0          "molecular weight of Kr")
(defparameter +M_Xe+  131.30d0         "molecular weight of Xe")
(defparameter +M_CH4+ 16.04303d0       "molecular weight of CH_4")
(defparameter +M_H2+  2.01594d0        "molecular weight of H_2")
(defparameter +N_A+   6.022169d26      "Avogadro's number [/kmole]")
(defparameter +Rstar+ 8.31432d3        "universal gas constant R^* [Nm/kmol.K]")

;; category II constants; km converted to m
(defparameter +F_N2+     0.78084d0     "fractional volume of N_2")
(defparameter +F_O2+     0.209476d0    "fractional volume of O_2")
(defparameter +F_Ar+     0.00934d0     "fractional volume of Ar")
(defparameter +F_CO2+    0.000314d0    "fractional volume of CO_2")
(defparameter +F_Ne+     0.00001818d0  "fractional volume of Ne")
(defparameter +F_He+     0.00000524d0  "fractional volume of He")
(defparameter +F_Kr+     0.00000114d0  "fractional volume of Kr")
(defparameter +F_Xe+     0.000000087d0 "fractional volume of Xe")
(defparameter +F_CH4+    0.000002d0    "fractional volume of CH_4")
(defparameter +F_H2+     0.0000005d0   "fractional volume of H_2")
(defparameter +g_0+      9.80665d0     "surface g [m/s^2]")
(defparameter +g_0prime+ +g_0+         "surface g\\prime [m^2/s.m\\prime]")
(defparameter +H_b+      #(0.0d3   11.0d3 20.0d3  32.0d3  47.0d3 51.0d3  71.0d3 84.8520d3)
  "reference levels [m]")
(defparameter +L_Mb+     #(-6.5d-3 0.0d-3 +1.0d-3 +2.8d-3 0.0d-3 -2.8d-3 -2.0d-3)
  "linear temperature gradients [K/m]")
(defparameter +P_0+      101325.0d0    "surface pressure [Pa]")
(defparameter +r_0+      6.356766d6    "radius of the Earth [m]")
(defparameter +T_0+      288.15d0      "surface temperature [K]")
(defparameter +S+        110.4d0       "Sutherland constant [K]; starndard p.4 gives 110, p.19 110.4")
(defparameter +beta+     1.458d-6      "[kg/s.m.K^{1/2}]; note error p4, correct p19")
(defparameter +gamma+    1.400d0       "ratio of specific heats")
(defparameter +sigma+    3.65d-10      "mean effective collision diameter [m]")

(defparameter +sigma+-sq (* +sigma+ +sigma+)) ;used mainly in this form

;; category III constants; km converted to m
;; because of Common Lisp case issue, lower case variables q_xxx, u_xxx and w_xxx
;; are replaced by qq_xxx, uu_xxx and ww_xxx
;; N_2 not present
(defparameter +a_O+      6.986d20) ;a_i in /m.s
(defparameter +a_O2+     4.863d20)
(defparameter +a_Ar+     4.487d20)
(defparameter +a_He+     1.700d21)
(defparameter +a_H+      3.305d21)
(defparameter +b_O+      0.750d0) ;b_i dimensionless
(defparameter +b_O2+     0.750d0)
(defparameter +b_Ar+     0.870d0)
(defparameter +b_He+     0.691d0)
(defparameter +b_H+      0.500d0)
(defparameter +K_7+      1.2d2) ;[m^2/s]
(defparameter +K_0+      0.0d0)
(defparameter +K_K7+     0.0d0) ;L_K,b in K/m, others 0
(defparameter +L_K9+     12.0d-3)
(defparameter +nO7+      8.6d16) ;[/m^3]
(defparameter +nH11+     8.0d10) ;[/m^3]
(defparameter +qq_O+     -3.416248d-12) ;/m^3; only for 86<=Z<97 km; for Z>97 km, +qq_O+ = 0
(defparameter +qq_O2+    0.0d0)
(defparameter +qq_Ar+    0.0d0)
(defparameter +qq_He+    0.0d0)
(defparameter +Q_O+      -5.809644d-13) ;[/m^3]
(defparameter +Q_O2+     1.366212d-13)
(defparameter +Q_Ar+     9.434079d-14)
(defparameter +Q_He+     -2.457369d-13)
(defparameter +T_9+      240.0d0)  ;[K]
(defparameter +T_inf+    1000.0d0) ;[K]
(defparameter +uu_O+     97.0d3)   ;[m]
(defparameter +uu_O2+    0.0d0)
(defparameter +uu_Ar+    0.0d0)
(defparameter +uu_He+    0.0d0)
(defparameter +U_O+      56.90311d3) ;[m]
(defparameter +U_O2+     86.000d3)
(defparameter +U_Ar+     86.000d3)
(defparameter +U_He+     86.000d3)
(defparameter +ww_O+     5.008765d-13) ;[/m^3]
(defparameter +ww_O2+    0.0d0)
(defparameter +ww_Ar+    0.0d0)
(defparameter +ww_He+    0.0d0)
(defparameter +W_O+      2.706240d-14) ;[/m^3]
(defparameter +W_O2+     8.333333d-14)
(defparameter +W_Ar+     8.333333d-14)
(defparameter +W_He+     6.666667d-13)
(defparameter +Z_b+      #(86.0d3 91.0d3 110.0d3 120.0d3 500.0d3 1000.0d3)) ;[m]; indexes are 7 to 12
(defparameter +alpha_N2+ 0.00d0) ;thermal diffusion coefficient for N_2; dimensionless
(defparameter +alpha_O+  0.00d0)
(defparameter +alpha_O2+ 0.00d0)
(defparameter +alpha_Ar+ 0.00d0)
(defparameter +alpha_He+ -0.40d0)
(defparameter +alpha_H+  -0.25d0)
(defparameter +phi+      7.2d11) ;[m^{-2}x^{-2} measure of vertical flux

;; calculated/derived constants
(defparameter +M_O+     (* 0.5d0 +M_O2+)
  "molecular weight of atomic oxygen") ;(letter O not zero)
(defparameter +M_H+     (* 0.5d0 +M_H2+)
  "molecular weight of atomic hydrogen")

(defparameter +M_0+ 28.9644d0 ;standard, section 1.2.4, p.9, gives this
  ;; calculation below gives 28.964507914535076 though
  ;; (/ (+ (* +F_N2+ +M_N2+) (* +F_O2+ +M_O2+) (* +F_Ar+ +M_Ar+)
  ;;       (* +F_CO2+ +M_CO2+) (* +F_Ne+ +M_Ne+) (* +F_He+ +M_He+)
  ;;       (* +F_Kr+ +M_Kr+) (* +F_Xe+ +M_Xe+) (* +F_CH4+ +M_CH4+) (* +F_H2+ +M_H2+))
  ;;    (+ +F_N2+ +F_O2+ +F_Ar+ +F_CO2+ +F_Ne+ +F_He+ +F_Kr+ +F_Xe+ +F_CH4+ +F_H2+))
  "average molecular weight of air at altitude zero")

(defparameter +T_Mb+
  (let ((T_Mb (make-array 8)))
    (setf (aref T_Mb 0) +T_0+)
    (loop
       for b from 1 to 7
       for b-1 = (- b 1) do
         (setf (aref T_Mb b)
               (+ (aref T_Mb b-1)
                  (* (aref +L_Mb+ b-1) (- (aref +H_b+ b) (aref +H_b+ b-1))))))
    T_Mb)
  "array of molecular-scale temperature at the bottom of each layer from 0 to 7")

(defparameter +T_10+
  (+ +T_9+ (* +L_K9+ (- (aref +Z_b+ 3) (aref +Z_b+ 2))))
  "kinetic temperature at Z_10")

(defparameter +lambda+
  (/ +L_K9+ (- +T_inf+ +T_10+))
  "used in formula for kinetic temperature between Z_10 to Z_12")

(defparameter +P_b+
  (let ((P_b (make-array 8)))
    (setf (aref P_b 0) +P_0+)
    (loop
       for b from 1 to 7
       for b-1 = (- b 1) do
         (setf (aref P_b b)
               (if (/= (aref +L_Mb+ b-1) 0.0d0)
                   (* (aref P_b b-1)
                      (expt (/ (aref +T_Mb+ b-1)
                               (+ (aref +T_Mb+ b-1)
                                  (* (aref +L_Mb+ b-1) (- (aref +H_b+ b) (aref +H_b+ b-1)))))
                            (/ (* +g_0prime+ +M_0+) (* +Rstar+ (aref +L_Mb+ b-1)))))
                   (* (aref P_b b-1)
                      (exp (/ (- (* +g_0prime+ +M_0+ (- (aref +H_b+ b) (aref +H_b+ b-1))))
                              (* +Rstar+ (aref +T_Mb+ b-1))))))))
    P_b)
  "array of pressure at the bottom of each layer from 0 to 7")


;;;; FULL ALTITUDE PROPERTIES

(defun Z-at-H (H)
  "Return geometric altitude corresponding to geopotential altitude, H, both in m."
  (/ (* +r_0+ H) (- +r_0+ H)))

(defun H-at-Z (Z)
  "Return geopotential altitude corresponding to geometric altitude, Z, both in m."
  (/ (* +r_0+ Z) (+ +r_0+ Z)))

(defun b-at-H (H)
  "Return zone number, b, 0..12, containing geopotential altitude H in m."
  (if (< H (aref +H_b+ 7))
      (loop for b from 0 to 6 do
           (when (< H (aref +H_b+ (+ b 1))) (return b))
         finally (return 6)) ;if return above didn't occur (should have)
      (loop with Z = (Z-at-H H)
         for b from 7 to 11 do
           (when (< Z (aref +Z_b+ (- b 6))) (return b))
         finally (return 12)))) ;or could issue a warning, beyond range

(defun b-at-Z (Z)
  "Return zone number, b, 0..12, containing geometric altitude, Z, in m."
  (if (>= Z (aref +Z_b+ 0))
      (loop for b from 7 to 11 do
           (when (< Z (aref +Z_b+ (- b 6))) (return b))
         finally (return 12)) ;or could issue a warning, beyond range
      (loop with H = (H-at-Z Z)
         for b from 0 to 6 do
           (when (< H (aref +H_b+ (+ b 1))) (return b))
         finally (return 6))))

(defun g-at-Z (Z)
  "Return the acceleration of gravity in m/s^2 at geometric altitude, Z, in m."
  (* +g_0+ (expt (/ +r_0+ (+ +r_0+ Z)) 2)))


;;;; Z <= 86 KM REGION PROPERTIES

;; Transition region from H = 79 km to H = 84.8520 km (Z = 86km) in which the molecular weight is
;; modified slightly to more smoothly transition between the lo and hi regions.
;; Standard, Table 8, p.9; [km] converted to [m].
(defparameter H-for-M-near-86km
  #(79.006d3 79.5d3 80.0d3 80.5d3
    81.0d3 81.5d3 82.0d3 82.5d3
    83.0d3 83.5d3 84.0d3 84.5d3 84.852d3))
(defparameter M/M_0-for-H-near-86km
  #(1.000000d0 0.999996d0 0.999988d0 0.999969
    0.999938d0 0.999904d0 0.999864d0 0.999822
    0.999778d0 0.999731d0 0.999681d0 0.999679d0 0.999579))

(defun M/M_0-at-H-near-86km (H)
  "Return M/M_0 for H values from 79.006 to 84.852 km (corresponding to Z values from 80 to 86 km).
Range is not checked."
  (interpolate H-for-M-near-86km M/M_0-for-H-near-86km H))

(defparameter Z-for-M-near-86km
  #(80.0d3 80.5d3 81.0d3 81.5d3
    82.0d3 82.5d3 83.0d3 83.5d3
    84.0d3 84.5d3 85.0d3 85.5d3 86.0d3))
(defparameter M/M_0-for-Z-near-86km
  #(1.000000d0 0.999996d0 0.999989d0 0.999971
    0.999941d0 0.999909d0 0.999870d0 0.999829
    0.999786d0 0.999741d0 0.999694d0 0.999641d0 0.999579))

(defun M/M_0-at-Z-near-86km (Z)
  "Return M/M_0 for Z values from 80 to 86 km.  Range is not checked."
  (interpolate-uniform 79.0d3 0.5d3 M/M_0-for-Z-near-86km Z))

(defun M/M_0-at-H-lo (H)
  "Return molecular weight at geopotential altitude H in m divided by the surface value of the
molecular weight.  Valid for H less than or equal to that value corresponding to geometric altitude 86
km (84.852d0 km)."
  (cond
    ((< H 79.006d3) 1.0d0)
    ((> H 84.852d3) 0.999579d0) ;catch values near 86 km that may be out due to
                        ;precision of 84.852d3 value
    (t (M/M_0-at-H-near-86km H))))

(defun M-at-H-lo (H)
  "Return molecular weight at geopotential altitude H in m.  Valid for H less than or equal to that
value corresponding to geometric altitude 86 km (84.852d0 km)."
  (cond
    ((< H 79.006d3) +M_0+)
    ((> H 84.852d3) (* +M_0+ 0.999579d0)) ;catch values near 86 km that may be out due to
                        ;precision of 84.852d3 value
    (t (* +M_0+ (M/M_0-at-H-near-86km H)))))

(defun T_M-at-H-lo-b (H b)
  "Return molecular-scale temperature at geopotential altitude H in m in zone B.
Does not check that H and B are consistent but assumes that they are.  Valid for H less than or equal
to that value corresponding to geometric altitude 86 km (84.852d0 km)."
  (let ((b (min b 6))) ;in  case of inconsitency around 86 km
    (+ (aref +T_Mb+ b) (* (aref +L_Mb+ b) (- H (aref +H_b+ b))))))

(defun T_M-at-H-lo (H)
  "Return molecular-scale temperature at geopotential altitude H in m.
Valid for H less than or equal to that value corresponding to geometric altitude 86 km (84.852d0 km)."
  (T_M-at-H-lo-b H (b-at-H H)))

(defun T-at-H-lo-b (H b)
  "Return kinetic temperature at geopotential altitude H in m in zone B.
Does not check that H and B are consistent but assumes that they are.  Valid for H less than or equal
to that value corresponding to geometric altitude 86 km (84.852d0 km)."
  (if (< H 79.006d3)
      (T_M-at-H-lo-b H b)
      (* (M/M_0-at-H-lo H) (T_M-at-H-lo-b H b))))

(defun T-at-H-lo (H)
  "Return kinetic temperature at geopotential altitude H in m.
Valid for H less than or equal to that value corresponding to geometric altitude 86 km (84.852d0 km)."
  (T-at-H-lo-b H (b-at-H H)))

(defun P-at-H-lo-b (H b)
  "Return pressure at geopotential altitude H in m in zone B.
Does not check that H and B are consistent but assumes that they are.  Valid for H less than or equal
to that value corresponding to geometric altitude 86 km (84.852d0 km)."
  (if (/= (aref +L_Mb+ b) 0.0d0)
      (* (aref +P_b+ b)
         (expt (/ (aref +T_Mb+ b) (+ (aref +T_Mb+ b) (* (aref +L_Mb+ b) (- H (aref +H_b+ b)))))
               (/ (* +g_0prime+ +M_0+) (* +Rstar+ (aref +L_Mb+ b)))))
      (* (aref +P_b+ b)
         (exp (/ (- (* +g_0prime+ +M_0+ (- H (aref +H_b+ b)))) (* +Rstar+ (aref +T_Mb+ b)))))))

(defun P-at-H-lo (H)
  "Return pressure at geopotential altitude H in m.  Valid for H less than or equal to that value
corresponding to geometric altitude 86 km (84.852d0 km)."
  (P-at-H-lo-b H (b-at-H H)))


;;;; Z > 86 KM REGION PROPERTIES

(defun T-at-Z-hi-b (Z b)
  "Return temperature in K at geometric altitude, Z, in meters in the high region, i.e., above 86 km.
Z is in zone B (not checked)."
  (case b
    (7 186.8673d0) ;this is with M/M_0 correction; 186.946d0 without
    (8 (let* ((T_c 263.1905d0)
              (A -76.3232d0)
              (aa -19.9429d3)
              (ZmZ8/a (/ (- Z 91.0d3) aa)))
         (+ T_c (* A (sqrt (- 1.0d0 (* ZmZ8/a ZmZ8/a)))))))
    (9 (+ +T_9+ (* +L_K9+ (- Z 110.0d3))))
    (t ; else for 10, 11, 12
     (let* ((Z_10 (aref +Z_b+ (- 10 7))) ;offset of 7 in +Z_b+ indexes
            (eps (/ (* (- Z Z_10) (+ +r_0+ Z_10)) (+ +r_0+ Z))))
       (- +T_inf+ (* (- +T_inf+ +T_10+) (exp (- (* +lambda+ eps)))))))))

(defun T-at-Z-hi (Z)
  "Return temperature in K at geometric altitude, Z, in meters in the high region, i.e., above 86 km."
  (T-at-Z-hi-b Z (b-at-Z Z)))

(defun dT/dZ-at-Z-b (Z b)
  "Return dT/dZ, the rate of change of kinetic temperature with altitude in K/m, at geometric altitude
Z in meters in the high region, i.e., above 86 km.  Z is in zone B (not checked)."
  (if (<= b 6)
      (aref +L_Mb+ b)
      (case b
        (7 0.0)
        (8 (let* ((A -76.3232)
                  (aa -19.9429d3) ;[m]
                  (ZmZ8/a (/ (- Z 91.0d3) aa)))
             (/ (* (/ (- A) aa) ZmZ8/a)
                (sqrt (- 1.0d0 (* ZmZ8/a ZmZ8/a))))))
        (9 12.0d-3)
        (t                               ;else for 10, 11, 12
         (let* ((Z_10 (aref +Z_b+ (- 10 7))) ;offset of 7 in +Z_b+ indexes
                (rZ (/ (+ +r_0+ Z_10) (+ +r_0+ Z)))
                (eps (* (- Z Z_10) rZ)))
           (* +lambda+ (- +T_inf+ +T_10+) rZ rZ (exp (- (* +lambda+ eps)))))))))

(defun dT/dZ-at-Z (Z)
  "Return dT/dZ, the rate of change of kinetic temperature with altitude in K/m,
at geometric altitude, Z, in meters in the high region, i.e., above 86 km."
  (dT/dZ-at-Z-b Z (b-at-Z Z)))

(defun f-integrand-N2 (Z g Tp)
  "Integrand for determining particle number distribution of N2 in the high region,
i.e., above 86 km."
  (let ((M (if (< Z 100.0d3) +M_0+ +M_N2+))) ;standard says <= (?)
    (/ (* M g) (* +Rstar+ Tp))))

(defun K (Z)
  "Return K value at Z.  See standard, eqns 7a, 7b, 7c."
  (cond ((< Z 95.0d3) +K_7+)
        ((< Z 115.0d3) (let ((Zm95km (- (* 1.0d-3 Z) 95.0d0)))
                         (* +K_7+ (exp (- 1.0d0 (/ 400.0 (- 400.0d0 (* Zm95km Zm95km))))))))
        (t 0.0d0)))

(defun f-integrand-O (Z g Tp n-background M-background)
  "Integrand for determining particle number distribution of O (atomic oxygen) in the high region,
i.e., above 86 km."
  ;; +alpha_O+ = 0
  (let* ((D_i (* (/ +a_O+ n-background) (expt (/ Tp 273.15d0) +b_O+))) ;[m^2/s]
         (K_Z (K Z))
         (fZ (if (= 0.0d0 K_Z)
                 (* (/ g (* +Rstar+ Tp)) (+ +M_O+))
                 (* (/ g (* +Rstar+ Tp)) (/ D_i (+ D_i K_Z))
                    (+ +M_O+ (* M-background (/ K_Z D_i))))))
         (ZmUi (- Z +U_O+))
         (ZmUi2 (* ZmUi ZmUi))
         (ZmUi3 (* ZmUi2 ZmUi))
         (uimZ (- +uu_O+ Z))
         (uimZ2 (* uimZ uimZ))
         (uimZ3 (* uimZ2 uimZ))
         (viOvDiK (+ (* +Q_O+ ZmUi2 (exp (- (* +W_O+ ZmUi3))))
                     (if (> Z 97000.0d0)
                         0.0
                         (* +qq_O+ uimZ2 (exp (- (* +ww_O+ uimZ3))))))))
    (+ fZ viOvDiK)))

(defun f-integrand-O2 (Z g Tp n-background M-background)
  "Integrand for determining particle number distribution of O2 in the high region,
i.e., above 86 km."
  ;; +alpha_O2+ = 0, +qq_O2+ = 0
  (let* ((D_i (* (/ +a_O2+ n-background) (expt (/ Tp 273.15d0) +b_O2+))) ;[m^2/s]
         (K_Z (K Z))
         (fZ (if (= 0.0d0 K_Z)
                 (* (/ g (* +Rstar+ Tp)) (+ +M_O2+))
                 (* (/ g (* +Rstar+ Tp)) (/ D_i (+ D_i K_Z))
                    (+ +M_O2+ (* M-background (/ K_Z D_i))))))
         (ZmUi (- Z +U_O2+))
         (ZmUi2 (* ZmUi ZmUi))
         (ZmUi3 (* ZmUi2 ZmUi))
         ;; (uimZ (- +uu_O2+ Z)) ;not used because q_i is 0
         ;; (uimZ2 (* uimZ uimZ))
         ;; (uimZ3 (* uimZ2 uimZ))
         (viOvDiK (* +Q_O2+ ZmUi2 (exp (- (* +W_O2+ ZmUi3))))))
    (+ fZ viOvDiK)))

(defun f-integrand-Ar (Z g Tp n-background M-background)
  "Integrand for determining particle number distribution of Ar in the high region,
i.e., above 86 km."
  ;; +alpha_Ar+ = 0, +qq_Ar+ = 0
  (let* ((D_i (* (/ +a_Ar+ n-background) (expt (/ Tp 273.15d0) +b_Ar+))) ;[m^2/s]
         (K_Z (K Z))
         (fZ (if (= 0.0d0 K_Z)
                 (* (/ g (* +Rstar+ Tp)) (+ +M_Ar+))
                 (* (/ g (* +Rstar+ Tp)) (/ D_i (+ D_i K_Z))
                    (+ +M_Ar+ (* M-background (/ K_Z D_i))))))
         (ZmUi (- Z +U_Ar+))
         (ZmUi2 (* ZmUi ZmUi))
         (ZmUi3 (* ZmUi2 ZmUi))
         ;; (uimZ (- +uu_Ar+ Z)) ;not used because q_i is 0
         ;; (uimZ2 (* uimZ uimZ))
         ;; (uimZ3 (* uimZ2 uimZ))
         (viOvDiK (* +Q_Ar+ ZmUi2 (exp (- (* +W_Ar+ ZmUi3))))))
    (+ fZ viOvDiK)))

(defun f-integrand-He (Z g Tp dT/dZ n-background M-background)
  "Integrand for determining particle number distribution of He in the high region,
i.e., above 86 km."
  ;; +qq_He+ = 0
  (let* ((D_i (* (/ +a_He+ n-background) (expt (/ Tp 273.15d0) +b_He+))) ;[m^2/s]
         (K_Z (K Z))
         (fZ (if (= 0.0d0 K_Z)
                 (* (/ g (* +Rstar+ Tp))
                    (+ +M_He+ (* (/ (* +alpha_He+ +Rstar+) g) dT/dZ)))
                 (* (/ g (* +Rstar+ Tp)) (/ D_i (+ D_i K_Z))
                    (+ +M_He+ (* M-background (/ K_Z D_i))
                       (* (/ (* +alpha_He+ +Rstar+) g) dT/dZ)))))
         (ZmUi (- Z +U_He+))
         (ZmUi2 (* ZmUi ZmUi))
         (ZmUi3 (* ZmUi2 ZmUi))
         ;; (uimZ (- +uu_He+ Z)) ;not used because q_i is 0
         ;; (uimZ2 (* uimZ uimZ))
         ;; (uimZ3 (* uimZ2 uimZ))
         (viOvDiK (* +Q_He+ ZmUi2 (exp (- (* +W_He+ ZmUi3))))))
    (+ fZ viOvDiK)))

(defparameter T_7
  186.8673d0
  "temperature at Z_7; this has M/M_0 correction which simple (aref +T_Mb+ 7) does not")

(defun make-v-params-for-n (Z-min Z-max dZ)
  "Make basic parameters used in particle number integration.  Return
 (values v-g v-Tp v-dT/dZ v-T7/T) for a uniform sequence of Z values
from Z-MIN to Z-MAX by DZ."
  (let* ((n (1+ (floor (- Z-max Z-min) dZ)))
         (v-g (make-array n))
         (v-Tp (make-array n))
         (v-dT/dZ (make-array n))
         (v-T7/T (make-array n)))
    (loop
       for i from 0 below n
       for Z = (+ Z-min (* i dZ))
       for Tp = (T-at-Z-hi Z) do
         (setf (aref v-g i) (g-at-Z Z))
         (setf (aref v-Tp i) Tp)
         (setf (aref v-dT/dZ i) (dT/dZ-at-Z Z))
         (setf (aref v-T7/T i) (/ T_7 Tp)))
    (values v-g v-Tp v-dT/dZ v-T7/T)))

(defun integrate-N2 (Z-start dZ v-g v-Tp v-T7/T)
  "Do N2 particle number integration for high region, i.e., above 86 km."
  (let* ((n_N2_7 1.129794d20)
         (n (length v-g))
         (v-integrand (make-array n))
         (v-n_N2 (make-array n)))
    (loop
       for i from 0 below n
       for Z = Z-start then (+ Z dZ) do
         (setf (aref v-integrand i) (f-integrand-N2 Z (aref v-g i) (aref v-Tp i))))
    (let ((v-integral (trapezoid-uniform dZ v-integrand)))
      (loop for i from 0 below n do
           (setf (aref v-n_N2 i) (* n_N2_7 (aref v-T7/T i) (exp (- (aref v-integral i)))))))
    v-n_N2))

(defun integrate-O-O2 (Z-start dZ v-g v-Tp v-T7/T v-n_N2)
  "Do O and O2 particle number integration for high region, i.e., above 86 km."
  (let* ((n_O_7 8.6d16)
         (n_O2_7 3.030898d19)
         (n (length v-g))
         (v-integrand-O (make-array n))
         (v-integrand-O2 (make-array n))
         (v-n_O (make-array n))
         (v-n_O2 (make-array n)))
    (loop
       for i from 0 below n
       for Z = Z-start then (+ Z dZ) do
         (setf (aref v-integrand-O i)
               (f-integrand-O Z (aref v-g i) (aref v-Tp i)
                              (aref v-n_N2 i) (if (< Z 100.0d3) +M_0+ +M_N2+)))
         (setf (aref v-integrand-O2 i)
               (f-integrand-O2 Z (aref v-g i) (aref v-Tp i)
                               (aref v-n_N2 i) (if (< Z 100.0d3) +M_0+ +M_N2+))))
    (let ((v-integral-O (trapezoid-uniform dZ v-integrand-O))
          (v-integral-O2 (trapezoid-uniform dZ v-integrand-O2)))
      (loop for i from 0 below n do
           (setf (aref v-n_O i) (* n_O_7 (aref v-T7/T i) (exp (- (aref v-integral-O i)))))
           (setf (aref v-n_O2 i) (* n_O2_7 (aref v-T7/T i) (exp (- (aref v-integral-O2 i)))))))
    (values v-n_O v-n_O2)))

(defun integrate-Ar-He (Z-start dZ v-g v-Tp v-dT/dZ v-T7/T v-n_N2 v-n_O v-n_O2)
  "Do Ar and He particle number integration for high region, i.e., above 86 km."
  (let* ((n_Ar_7 1.351400d18)
         (n_He_7 7.5817d14)
         (n (length v-g))
         (v-integrand-Ar (make-array n))
         (v-integrand-He (make-array n))
         (v-n_Ar (make-array n))
         (v-n_He (make-array n)))
    (loop
       for i from 0 below n
       for Z = Z-start then (+ Z dZ) do
         (let* ((n_N2 (aref v-n_N2 i))
                (n_O (aref v-n_O i))
                (n_O2 (aref v-n_O2 i))
                (n-bkg (+ n_N2 n_O n_O2))
                (M-bkg (if (<= Z 100.0d3)
                           +M_0+
                           (/ (+ (* n_N2 +M_N2+) (* n_O +M_O+) (* n_O2 +M_O2+))
                              (+ n_N2 n_O n_O2)))))
           (setf (aref v-integrand-Ar i)
                 (f-integrand-Ar Z (aref v-g i) (aref v-Tp i) n-bkg M-bkg))
           (setf (aref v-integrand-He i)
                 (f-integrand-He Z (aref v-g i) (aref v-Tp i) (aref v-dT/dZ i) n-bkg M-bkg))))
    (let ((v-integral-Ar (trapezoid-uniform dZ v-integrand-Ar))
          (v-integral-He (trapezoid-uniform dZ v-integrand-He)))
      (loop for i from 0 below n do
           (setf (aref v-n_Ar i) (* n_Ar_7 (aref v-T7/T i) (exp (- (aref v-integral-Ar i)))))
           (setf (aref v-n_He i) (* n_He_7 (aref v-T7/T i) (exp (- (aref v-integral-He i)))))))
    (values v-n_Ar v-n_He)))

(defun integrate-H (Z-min Z-tau Z11 Z-max dZ
                    v-g v-Tp v-n_N2 v-n_O v-n_O2 v-n_Ar v-n_He)
  "Do H (atomic hydrogen) particle number integration for high region, i.e., above 86 km.
Assume uniform Z spacing in input arrays.
Z-MIN is the height at which the input arrays begin, and the output arrays begin.
Z-TAU is the height above which atomic hydrogen is calculated (150 km in the spec).
It is zero from Z-MIN to below Z-TAU."
  (let* ((T11 999.2356d0)
         (alphaH+1 (+ 1.0d0 -0.25d0))
         (i-at-Z-min 0)
         (i-at-Z-tau (floor (- Z-tau Z-min) dZ)) ;i-tau is 0 at Z-tau where i = i-at-Z-tau
         (i-at-Z11 (floor (- Z11 Z-min) dZ))
         (i-at-Z-max (1+ (floor (- Z-max Z-min) dZ))) ;allow one extra in case of roundoff
         (n-Z11-down (1+ (- i-at-Z11 i-at-Z-tau))) ;1+ to include Z11
         (n-Z11-up (- i-at-Z-max i-at-Z11)) ;does not include Z11
         (v-Z-Z11-down (make-array n-Z11-down))
         (v-Z-Z11-up (make-array n-Z11-up))
         (v-n_H (make-array (1+ (- i-at-Z-max i-at-Z-min))))
         (M_H/Rstar (/ +M_H+ +Rstar+)))
    (loop for i from (1- n-Z11-down) downto 0
       for Z = Z11 then (- Z dZ)
       do (setf (aref v-Z-Z11-down i) Z))
    (loop for i from 0 below n-Z11-up
       for Z = Z11 then (+ Z dZ)
       do (setf (aref v-Z-Z11-up i) Z))
    (let ((v-tau-Z11-down)
          (v-tau-Z11-up)
          ;;(v-n_H-Z11-down)
          (v-integral-H-Z11-down))
      ;; tau integration from Z11 down
      (let* ((v-integrand (make-array n-Z11-down)))
        (loop for i-tau from 0 below n-Z11-down
           for i = (+ i-at-Z-tau i-tau)
           do (setf (aref v-integrand i-tau) (/ (* (aref v-g i) M_H/Rstar) (aref v-Tp i))))
        (setf v-tau-Z11-down (trapezoid-uniform dZ v-integrand :start (1- n-Z11-down) :end 0)))
      ;; tau integration from Z11 up
      (let* ((v-integrand (make-array n-Z11-up)))
        (loop for i-tau from 0 below n-Z11-up
           for i = (+ i-at-Z11 i-tau) do
             (setf (aref v-integrand i-tau) (/ (* (aref v-g i) M_H/Rstar) (aref v-Tp i))))
        (setf v-tau-Z11-up (trapezoid-uniform dZ v-integrand)))
      ;; H integration Z11 down
      (let* ((v-integrand (make-array n-Z11-down)))
        (loop for i-tau from 0 below n-Z11-down
           for i = (+ i-at-Z-tau i-tau)
           for Tp = (aref v-Tp i)
           for n-bkg = (+ (aref v-n_N2 i) (aref v-n_O i) (aref v-n_O2 i)
                          (aref v-n_Ar i) (aref v-n_He i))
           for D_H = (* (/ +a_H+ n-bkg) (expt (/ Tp 273.15d0) +b_H+)) do
             (setf (aref v-integrand i-tau)
                   (* (/ +phi+ D_H) (expt (/ Tp T11) alphaH+1) (exp (aref v-tau-Z11-down i-tau)))))
        (setf v-integral-H-Z11-down (trapezoid-uniform dZ v-integrand :start (1- n-Z11-down) :end 0)))
      ;; make final H
      (loop for i from 0 below i-at-Z-tau do
           (setf (aref v-n_H i) 0.0d0))
      (loop for i-tau from 0 below n-Z11-down
         for i = (+ i-at-Z-tau i-tau)
         for Tp = (aref v-Tp i) do
           (setf (aref v-n_H i) (* (- +nH11+ (aref v-integral-H-Z11-down i-tau))
                                   (expt (/ T11 Tp) alphaH+1)
                                   (exp (- (aref v-tau-Z11-down i-tau))))))
      (loop             ;H Z11 up in diffusive equilibrium
         for i-tau from 0 below n-Z11-up
         for i = (+ i-at-Z11 i-tau)
         for Tp = (aref v-Tp i) do
           (setf (aref v-n_H i) (* (* +nH11+ (expt (/ T11 Tp) alphaH+1))
                                   (exp (- (aref v-tau-Z11-up i-tau))))))
      v-n_H)))


;;;; USSADATA STRUCTURE

(defstruct ussadata
  ;; structure containing US Standard Atmosphere data for the region above 86 km for later
  ;; interpolation
  Z-start-hi ;starting Z for uniformly spaced Z values corresponding to the vectors below
  dZ-hi      ;Z spacing for the vectors below
  v-n_N2-hi  ;vector of N2 particle number density
  v-n_O-hi   ;vector of O (atomic oxygen) particle number density
  v-n_O2-hi  ;vector of O2 particle number density
  v-n_Ar-hi  ;vector of Ar particle number density
  v-n_He-hi  ;vector of He particle number density
  v-n_H-hi   ;vector of H (atomic hydrogen) particle number density
  v-n-hi     ;vector of total number density
  v-sumMn-hi ;vector of sum-over-species of M*n
  v-M-hi     ;vector of molecular weights
  )

(defun make-ussadata-with-dZ-hi (&optional (dZ-hi 10.0d0))
                        ;dZ-hi = 10 m gives close match to the standard but is probably much
                        ;smaller than would be needed for realistic tolerance
  "Make a ussadata structure holding US Standard Atmosphere data for later interpolation.  Data in the
high region, above 86 km, is calculated and stored for altitudes separated by DZ-HI."
  (let ((Z-min 86.0d3)
        (Z-tau 150.0d3)
        (Z11 500.0d3)
        (Z-max 1000.0d3))
    (multiple-value-bind (v-g v-Tp v-dT/dZ v-T7/T) (make-v-params-for-n Z-min Z-max dZ-hi)
      (let* ((v-n_N2 (integrate-N2 Z-min dZ-hi v-g v-Tp v-T7/T)))
        (multiple-value-bind (v-n_O v-n_O2) (integrate-O-O2 Z-min dZ-hi v-g v-Tp v-T7/T v-n_N2)
          (multiple-value-bind (v-n_Ar v-n_He)
              (integrate-Ar-He Z-min dZ-hi v-g v-Tp v-dT/dZ v-T7/T v-n_N2 v-n_O v-n_O2)
            (let* ((v-n_H (integrate-H Z-min Z-tau Z11 Z-max dZ-hi
                                       v-g v-Tp v-n_N2 v-n_O v-n_O2 v-n_Ar v-n_He))
                   (v-n (map 'vector #'+ v-n_N2 v-n_O v-n_O2 v-n_Ar v-n_He v-n_H))
                   (v-sumMn (map 'vector (lambda (n_N2 n_O n_O2 n_Ar n_He n_H)
                                           (+ (* n_N2 +M_N2+) (* n_O +M_O+) (* n_O2 +M_O2+)
                                              (* n_Ar +M_Ar+) (* n_He +M_He+) (* n_H +M_H+)))
                                 v-n_N2 v-n_O v-n_O2 v-n_Ar v-n_He v-n_H))
                   (v-M (map 'vector #'/ v-sumMn v-n)))
              (make-ussadata :Z-start-hi Z-min :dZ-hi dZ-hi
                             :v-n_N2-hi v-n_N2 :v-n_O-hi v-n_O :v-n_O2-hi v-n_O2
                             :v-n_Ar-hi v-n_Ar :v-n_He-hi v-n_He :v-n_H-hi v-n_H
                             :v-n-hi v-n
                             :v-sumMn-hi v-sumMn :v-M-hi v-M))))))))

(defparameter *ussadata* (make-ussadata-with-dZ-hi)
  "Default ussadata data object.  API functions interpolate the data stored in this object in the high
region, i.e., above 86 km.")

(defmacro with-ussadata (ussadata &rest body)
  "Execute the forms in BODY with *ussadata* bound to USSADATA.  Restore *ussadata* to its
original value when this is done."
  (let ((current-ussadata (gensym)))
    `(let ((,current-ussadata *ussadata*))
       (setf *ussadata* ,ussadata)
       (let ((ret (progn ,@body)))
         (setf *ussadata* ,current-ussadata)
         ret))))


;;;; PROPERTIES AT GEOMETRIC ALTITUDE Z

;; H-at-Z is above

(defun M/M_0-at-Z (Z)
  "Return M/M_0, the ratio of molecular weight at geometric altitude, Z, in meters,
to the molecular weight at altitude zero."
  (cond ((< Z 80.0d3) 1.0d0)
        ((< Z 86.0d3) (M/M_0-at-Z-near-86km Z))
        (t (/ (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                   (ussadata-v-M-hi *ussadata*) Z) +M_0+))))

(defun M-at-Z (Z)
  "Return M/M_0, the molecular weight at geometric altitude, Z, in meters."
  (cond ((< Z 80.0d3) +M_0+)
        ((< Z 86.0d3) (* +M_0+ (M/M_0-at-Z-near-86km Z)))
        (t (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                (ussadata-v-M-hi *ussadata*) Z))))

(defun T-at-Z (Z)
  "Return the kinetic temperature in K at geometric altitude, Z, in meters."
  (cond ((< Z 86.0d3) (T-at-H-lo (H-at-Z Z)))
        (t (T-at-Z-hi Z))))

(defun T_M-at-Z (Z)
  "Return the molecular-scale temperature in K at geometric altitude, Z, in meters."
  (cond ((< Z 86.0d3) (T_M-at-H-lo (H-at-Z Z)))
        (t (/ (T-at-Z Z) (M/M_0-at-Z Z))) ;strictly speaking, not defined
        ))

(defun n-at-Z (Z)
  "Return the total particle number density in /m^3 at geometric altitude, Z, in meters."
  (cond ((< Z 80.0d3) (let* ((H (H-at-Z Z))
                             (b (b-at-H H))
                             (P (P-at-H-lo-b H b))
                             (Tp (T_M-at-H-lo-b H b))) ;T_M = T in this region
                        (/ (* +N_A+ P) (* Tp +Rstar+))))
        ((< Z 86.0d3) (let* ((H (H-at-Z Z))
                             (b (b-at-H H))
                             (P (P-at-H-lo-b H b))
                             (T_M (T-at-H-lo-b H b))
                             (M/M_0 (M/M_0-at-H-lo H)))
                        (/ (* +N_A+ P) (* T_M +Rstar+ M/M_0))))
        (t (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                (ussadata-v-n-hi *ussadata*) Z))))

(defun n_N2-at-Z (Z)
  "Return number density of N2 in /m^3 at geometric altitude, Z, in meters."
  (cond ((< Z 86.0d3)
         (let ((n (n-at-Z Z))) (* n +F_N2+)))
        (t (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                (ussadata-v-n_N2-hi *ussadata*) Z))))

(defun n_O-at-Z (Z)
  "Return number density of O in /m^3 at geometric altitude, Z, in meters."
  (cond ((< Z 86.0d3) 0.0d0)
        (t (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                (ussadata-v-n_O-hi *ussadata*) Z))))

(defun n_O2-at-Z (Z)
  "Return number density of O2 in /m^3 at geometric altitude, Z, in meters."
  (cond ((< Z 86.0d3)
         (let ((n (n-at-Z Z))) (* n +F_O2+)))
        (t (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                (ussadata-v-n_O2-hi *ussadata*) Z))))

(defun n_Ar-at-Z (Z)
  "Return number density of Ar in /m^3 at geometric altitude, Z, in meters."
  (cond ((< Z 86.0d3)
         (let ((n (n-at-Z Z))) (* n +F_Ar+)))
        (t (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                (ussadata-v-n_Ar-hi *ussadata*) Z))))

(defun n_He-at-Z (Z)
  "Return number density of He in /m^3 at geometric altitude, Z, in meters."
  (cond ((< Z 86.0d3)
         (let ((n (n-at-Z Z))) (* n +F_He+)))
        (t (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                (ussadata-v-n_He-hi *ussadata*) Z))))

(defun n_H-at-Z (Z)
  "Return number density of H in /m^3 at geometric altitude, Z, in meters."
  (cond ((< Z 86.0d3) 0.0d0)
        (t (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                (ussadata-v-n_H-hi *ussadata*) Z))))

(defun v-n_X-at-Z (Z)
  "Return a 7-vector of number densities for N2, O, O2, Ar, He, H, and the total
in /m^3 at geometric altitude, Z, in meters."
  (cond ((< Z 86.0d3) (let ((n (n-at-Z Z)))
                        (vector (* n +F_N2+) 0.0d0 (* n +F_O2+) (* n +F_Ar+) (* n +F_He+) 0.0d0 n)))
        (t (multiple-value-bind (span a1 a2)
               (interpolate-uniform-weights
                (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                (ussadata-v-n_N2-hi *ussadata*) Z)
             (let* ((v-n_N2 (ussadata-v-n_N2-hi *ussadata*))
                    (v-n_O (ussadata-v-n_O-hi *ussadata*))
                    (v-n_O2 (ussadata-v-n_O2-hi *ussadata*))
                    (v-n_Ar (ussadata-v-n_Ar-hi *ussadata*))
                    (v-n_He (ussadata-v-n_He-hi *ussadata*))
                    (v-n_H (ussadata-v-n_H-hi *ussadata*))
                    (v-n (ussadata-v-n-hi *ussadata*))
                    (n_N2 (+ (* a1 (aref v-n_N2 span)) (* a2 (aref v-n_N2 (1+ span)))))
                    (n_O (+ (* a1 (aref v-n_O span)) (* a2 (aref v-n_O (1+ span)))))
                    (n_O2 (+ (* a1 (aref v-n_O2 span)) (* a2 (aref v-n_O2 (1+ span)))))
                    (n_Ar (+ (* a1 (aref v-n_Ar span)) (* a2 (aref v-n_Ar (1+ span)))))
                    (n_He (+ (* a1 (aref v-n_He span)) (* a2 (aref v-n_He (1+ span)))))
                    (n_H (+ (* a1 (aref v-n_H span)) (* a2 (aref v-n_H (1+ span)))))
                    (n (+ (* a1 (aref v-n span)) (* a2 (aref v-n (1+ span))))))
               (vector n_N2 n_O n_O2 n_Ar n_He n_H n))))))

(defun P-at-Z (Z)
  "Return pressure in Pa at geometric altitude, Z, in meters"
  (cond ((< Z 86.0d3) (P-at-H-lo (H-at-Z Z)))
        (t (* (/ (* +Rstar+ (T-at-Z-hi Z)) +N_A+)
              (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                   (ussadata-v-n-hi *ussadata*) Z)))))

(defparameter +P_0+ (P-at-Z 0.0d0) "Pressure at altitude 0.")

(defun P/P_0-at-Z (Z)
  "Return pressure at geometric altitude, Z, in meters divided by the pressure at altitude zero."
  (/ (P-at-Z Z) +P_0+))

(defun rho-at-Z (Z)
  "Return density in kg/m^3 at geometric altitude, Z, in meters."
  (cond ((< Z 80.0d3) (let* ((H (H-at-Z Z))
                             (b (b-at-H H))
                             (P (P-at-H-lo-b H b))
                             (T_M (T_M-at-H-lo-b H b)))
                        (/ (* P +M_0+) (* +Rstar+ T_M))))
        ((< Z 86.0d3) (let* ((H (H-at-Z Z))
                             (b (b-at-H H))
                             (P (P-at-H-lo-b H b))
                             (T_M (T_M-at-H-lo-b H b)))
                        (/ (* P +M_0+) (* +Rstar+ T_M)))) ;same actually
        (t (/ (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                   (ussadata-v-sumMn-hi *ussadata*) Z) +N_A+))))

(defparameter +rho_0+ (rho-at-Z 0.0d0)
  "Density at altitude 0.")

(defun rho/rho_0-at-Z (Z)
  "Return density at geometric altitude, Z, in meters divided by the density at altitude zero."
  (/ (rho-at-Z Z) +rho_0+))

(defun v_m-at-Z (Z)
  "Return mole volume in m^3 at geometric altitude, Z, in meters."
  (cond ((< Z 80.0d3) (let* ((H (H-at-Z Z))
                             (b (b-at-H H))
                             (P (P-at-H-lo-b H b))
                             (Tp (T_M-at-H-lo-b H b))) ;T = T_M in this region
                        (/ (* +Rstar+ Tp) P)))
        ((< Z 86.0d3) (let* ((H (H-at-Z Z))
                             (b (b-at-H H))
                             (P (P-at-H-lo-b H b))
                             (Tp (T-at-H-lo-b H b))) ;T != T_M in this region
                        (/ (* +Rstar+ Tp) P)))
        (t (/ +N_A+ (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                         (ussadata-v-n-hi *ussadata*) Z)))))

(defun H_P-at-Z (Z)
  "Return pressure scale height in meters at geometric altitude, Z, in meters."
  (cond ((< Z 86.0d3) (let* ((H (H-at-Z Z))
                             (b (b-at-H H))
                             (T_M (T_M-at-H-lo-b H b))
                             (g (g-at-Z Z)))
                        (/ (* +Rstar+ T_M) (* g +M_0+))))
        (t (let* ((Tp (T-at-Z-hi Z))
                  (n (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                          (ussadata-v-n-hi *ussadata*) Z))
                  (sumMn (interpolate-uniform
                          (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                          (ussadata-v-sumMn-hi *ussadata*) Z))
                  (g (g-at-Z Z)))
             ;; does not quite match the standard, though its individual components do (?);
             ;; for example at 10000 m, T_M = T = 223.25209264797857d0, g = 9.775868442887434d0,
             ;; and the constants +Rstar+ = 8314.32d0 and +M_0+ = 28.9644d0;
             ;; this gives H_P = 6555.448184469401d0 but the standard gives 6552.5 ??
             (/ (* +Rstar+ Tp n) (* g sumMn))))))

(defun H_rho-at-Z (Z &key (dZ 100.0d0))
  "Return density scale height in meters at geometric altitude, Z, in meters.  In the high altitude region,
above 86 km, this uses a numerical derivative to get d ln(T_M) / d Z.  DZ is the half-height over
which this derivative is done, centered around Z."
  ;; may as well just do numerical differentiation of ln(rho) then !
  (cond ((< Z 86.0d3) (let* ((H_P (H_P-at-Z Z))
                             (dT/dZ (dT/dZ-at-Z Z))
                             (Tp (T-at-Z Z))
                             (dlnTdZ (/ dT/dZ Tp)))
                        (/ H_P (+ 1.0d0 (* H_P dlnTdZ)))))
        (t (let* ((H_P (H_P-at-Z Z))
                  (T_M- (T_M-at-Z (- Z dZ)))
                  (T_M+ (T_M-at-Z (+ Z dZ)))
                  (dlnT_MdZ (/ (- (log T_M+) (log T_M-)) (* 2.0d0 dZ))))
             (/ H_P (+ 1.0d0 (* H_P dlnT_MdZ)))))))

(defun V-at-Z (Z)
  "Return mean air particle speed in m/s at geometric altitude, Z, in meters."
  (cond ((< Z 86.0d3) (let* ((H (H-at-Z Z))
                             (T_M (T_M-at-H-lo H)))
                        (sqrt (/ (* 8.0d0 +Rstar+ T_M) pi +M_0+))))
        (t (let ((Tp (T-at-Z-hi Z))
                 (n (interpolate-uniform
                     (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                     (ussadata-v-n-hi *ussadata*) Z))
                 (sumMn (interpolate-uniform
                         (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                         (ussadata-v-sumMn-hi *ussadata*) Z)))
             (sqrt (/ (* 8.0d0 +Rstar+ Tp n) pi sumMn))))))

(defun L-at-Z (Z)
  "Return mean free path in meters at geometric altitude, Z, in meters."
  (cond ((< Z 80.0d3) (let* ((H (H-at-Z Z))
                             (b (b-at-H H))
                             (T_M (T_M-at-H-lo-b H b))
                             (P (P-at-H-lo-b H b)))
                        (/ (* +Rstar+ T_M) (* (sqrt 2.0d0) pi +N_A+ +sigma+-sq P))))
        ((< Z 86.0d3) (let* ((H (H-at-Z Z))
                             (b (b-at-H H))
                             (T_M (T_M-at-H-lo-b H b))
                             (M/M_0 (M/M_0-at-Z-near-86km Z))
                             (P (P-at-H-lo-b H b)))
                        (/ (* +Rstar+ M/M_0 T_M)
                           (* (sqrt 2.0d0) pi +N_A+ +sigma+-sq P))))
        (t (let ((n (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                         (ussadata-v-n-hi *ussadata*) Z)))
             (/ (* (sqrt 2.0d0) pi +sigma+-sq n))))))

(defun nu-at-Z (Z)
  "Return mean collision frequency in /s at geometric altitude, Z, in meters."
  (cond ((< Z 80.0d3) (let* ((H (H-at-Z Z))
                             (b (b-at-H H))
                             (T_M (T_M-at-H-lo-b H b))
                             (P (P-at-H-lo-b H b)))
                        (* 4.0d0 +N_A+ +sigma+-sq
                           (sqrt (/ (* pi P P) (* +Rstar+ +M_0+ T_M))))))
        ((< Z 86.0d3) (let* ((H (H-at-Z Z))
                             (b (b-at-H H))
                             (T_M (T_M-at-H-lo-b H b))
                             (M/M_0 (M/M_0-at-Z-near-86km Z))
                             (M (* +M_0+ M/M_0))
                             (P (P-at-H-lo-b H b)))
                        (* 4.0d0 +N_A+ +sigma+-sq
                           (sqrt (/ (* pi P P) (* +Rstar+ M T_M))))))
        (t (let* ((n (interpolate-uniform
                      (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                      (ussadata-v-n-hi *ussadata*) Z))
                  (sumMn (interpolate-uniform
                          (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                          (ussadata-v-sumMn-hi *ussadata*) Z))
                  (P (P-at-Z Z))
                  (Tp (T-at-Z-hi Z)))
             (* 4.0d0 +N_A+ +sigma+-sq
                (sqrt (/ (* pi P P n) (* +Rstar+ Tp sumMn))))))))

(defun C_s-at-Z (Z)
  "Return speed of sound in m/s at geometric altitude, Z, in meters."
  (sqrt (/ (* +gamma+ +Rstar+ (T_M-at-Z Z)) +M_0+)))

(defun mu-at-Z (Z)
  "Return dynamic viscosity in N.s/m^2 at geometric altitude, Z, in meters."
  (let ((Tp (T-at-Z Z)))
    (/ (* +beta+ (expt Tp 1.5d0)) (+ Tp +S+))))

(defparameter +mu_0+ (mu-at-Z 0.0d0)
  "Dynamic viscosity at latitude zero.")

(defun mu/mu_0-at-Z (Z)
  "Return dynamic viscosity at geometric altitude, Z, in meters divided by the
dynamic viscosity at altitude zero."
  (/ (mu-at-Z Z) +mu_0+))

(defun eta-at-Z (Z)
  "Return kinematic viscosity in m^2/s at geometric altitude, Z, in meters."
  (/ (mu-at-Z Z) (rho-at-Z Z)))

(defparameter +eta_0+ (eta-at-Z 0.0d0)
  "Kinematic viscosity at altitude zero.")

(defun eta/eta_0-at-Z (Z)
  "Return kinematic viscosity at geometric altitude, Z, in meters divided by the
kinematic viscosity at altitude zero."
  (/ (eta-at-Z Z) +eta_0+))

(defun k_t-at-Z (Z)
  "Return thermal conductivity in W/m.K at geometric altitude, Z, in meters."
  (let ((Tp (T-at-Z Z)))
    (/ (* 2.65019d-3 (expt Tp 1.5d0))
       (+ Tp (* 245.4d0 (expt 10.0d0 (/ -12.0d0 Tp)))))))

(defparameter +k_t0+ (k_t-at-Z 0.0d0)
  "Thermal conductivity at altitude zero.")

(defun k_t/k_t0-at-Z (Z)
  "Return thermal conductivity at geometric altitude, Z, in meters divided by the
thermal conductivity at altitude zero."
  (/ (k_t-at-Z Z) +k_t0+))


;;;; PROPERTIES AT GEOPOTENTIAL ALTITUDE H

;; Z-at-H is above

(defun g-at-H (H)
  "Return the acceleration of gravity in m/s^2 at geopotential altitude, H, in m."
  (g-at-Z (Z-at-H H)))

(defun M/M_0-at-H (H)
  "Return M/M_0, the ratio of molecular weight at geopotential altitude, H, in meters,
to the molecular weight at altitude zero."
  (cond ((< H 84.852d3) (M/M_0-at-H-lo H))
        (t (/ (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                   (ussadata-v-M-hi *ussadata*) (Z-at-H H))
              +M_0+))))

(defun M-at-H (H)
  "Return M/M_0, the molecular weight at geopotential altitude, H, in meters."
  (cond ((< H 84.852d3) (M-at-H-lo H))
        (t (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                (ussadata-v-M-hi *ussadata*) (Z-at-H H)))))

(defun T-at-H (H)
  "Return the kinetic temperature in K at geopotential altitude, H, in meters."
  (cond ((< H 84.852d3) (T-at-H-lo H))
        (t (T-at-Z-hi (Z-at-H H)))))

(defun T_M-at-H (H)
  "Return the molecular-scale temperature in K at geopotential altitude, H, in meters."
  (cond ((< H 84.852d3) (T_M-at-H-lo H))
        (t (let ((Z (Z-at-H H))) (/ (T-at-Z Z) (M/M_0-at-Z Z)))) ;strictly speaking, not defined
        ))

(defun n-at-H (H)
  "Return the total particle number density in /m^3 at geopotential altitude, H, in meters."
  (cond ((< H 79.006d3) (let* ((b (b-at-H H))
                               (P (P-at-H-lo-b H b))
                               (Tp (T_M-at-H-lo-b H b))) ;T_M = T in this region
                          (/ (* +N_A+ P) (* Tp +Rstar+))))
        ((< H 84.852d3) (let* ((b (b-at-H H))
                               (P (P-at-H-lo-b H b))
                               (T_M (T-at-H-lo-b H b))
                               (M/M_0 (M/M_0-at-H-lo H)))
                          (/ (* +N_A+ P) (* T_M +Rstar+ M/M_0))))
        (t (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                (ussadata-v-n-hi *ussadata*) (Z-at-H H)))))

(defun n_N2-at-H (H)
  "Return number density of N2 in /m^3 at geopotential altitude, H, in meters."
  (cond ((< H 84.852d3) (* (n-at-H H) +F_N2+))
        (t (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                (ussadata-v-n_N2-hi *ussadata*) (Z-at-H H)))))

(defun n_O-at-H (H)
  "Return number density of O in /m^3 at geopotential altitude, H, in meters."
  (cond ((< H 84.852d3) 0.0d0)
        (t (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                (ussadata-v-n_O-hi *ussadata*) (Z-at-H H)))))

(defun n_O2-at-H (H)
  "Return number density of O2 in /m^3 at geopotential altitude, H, in meters."
  (cond ((< H 84.852d3) (* (n-at-H H) +F_O2+))
        (t (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                (ussadata-v-n_O2-hi *ussadata*) (Z-at-H H)))))

(defun n_Ar-at-H (H)
  "Return number density of Ar in /m^3 at geopotential altitude, H, in meters."
  (cond ((< H 84.852d3) (* (n-at-H H) +F_Ar+))
        (t (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                (ussadata-v-n_Ar-hi *ussadata*) (Z-at-H H)))))

(defun n_He-at-H (H)
  "Return number density of He in /m^3 at geopotential altitude, H, in meters."
  (cond ((< H 84.852d3) (* (n-at-H H) +F_He+))
        (t (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                (ussadata-v-n_He-hi *ussadata*) (Z-at-H H)))))

(defun n_H-at-H (H)
  "Return number density of H in /m^3 at geopotential altitude, H, in meters."
  (cond ((< H 84.852d3) 0.0d0)
        (t (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                (ussadata-v-n_H-hi *ussadata*) (Z-at-H H)))))

(defun v-n_X-at-H (H)
  "Return a 7-vector of number densities for N2, O, O2, Ar, He, H, and the total
in /m^3 at geopotential altitude, H, in meters."
  (cond ((< H 84.852d3) (let ((n (n-at-Z (Z-at-H H))))
                          (vector (* n +F_N2+) 0.0d0 (* n +F_O2+) (* n +F_Ar+) (* n +F_He+) 0.0d0 n)))
        (t (multiple-value-bind (span a1 a2)
               (interpolate-uniform-weights
                (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                            (ussadata-v-n_N2-hi *ussadata*) (Z-at-H H))
             (let* ((v-n_N2 (ussadata-v-n_N2-hi *ussadata*))
                    (v-n_O (ussadata-v-n_O-hi *ussadata*))
                    (v-n_O2 (ussadata-v-n_O2-hi *ussadata*))
                    (v-n_Ar (ussadata-v-n_Ar-hi *ussadata*))
                    (v-n_He (ussadata-v-n_He-hi *ussadata*))
                    (v-n_H (ussadata-v-n_H-hi *ussadata*))
                    (v-n (ussadata-v-n-hi *ussadata*))
                    (n_N2 (+ (* a1 (aref v-n_N2 span)) (* a2 (aref v-n_N2 (1+ span)))))
                    (n_O (+ (* a1 (aref v-n_O span)) (* a2 (aref v-n_O (1+ span)))))
                    (n_O2 (+ (* a1 (aref v-n_O2 span)) (* a2 (aref v-n_O2 (1+ span)))))
                    (n_Ar (+ (* a1 (aref v-n_Ar span)) (* a2 (aref v-n_Ar (1+ span)))))
                    (n_He (+ (* a1 (aref v-n_He span)) (* a2 (aref v-n_He (1+ span)))))
                    (n_H (+ (* a1 (aref v-n_H span)) (* a2 (aref v-n_H (1+ span)))))
                    (n (+ (* a1 (aref v-n span)) (* a2 (aref v-n (1+ span))))))
               (vector n_N2 n_O n_O2 n_Ar n_He n_H n)))
           ;; (interpolate-v (ussadata-v-Z-hi *ussadata*) (ussadata-v-v-n-hi *ussadata*) (Z-at-H H))
           )))

(defun P-at-H (H)
  "Return pressure in Pa at geopotential altitude, H, in meters"
  (cond ((< H 84.852d3) (P-at-H-lo H))
        (t (let ((Z (Z-at-H H)))
             (* (/ (* +Rstar+ (T-at-Z-hi Z)) +N_A+)
                (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                     (ussadata-v-n-hi *ussadata*) Z))))))

(defun P/P_0-at-H (H)
  "Return pressure at geopotential altitude, H, in meters divided by the pressure at altitude zero."
  (/ (P-at-H H) +P_0+))

(defun rho-at-H (H)
  "Return density in kg/m^3 at geopotential altitude, H, in meters."
  (cond ((< H 79.006d3) (let* ((b (b-at-H H))
                               (P (P-at-H-lo-b H b))
                               (T_M (T_M-at-H-lo-b H b)))
                          (/ (* P +M_0+) (* +Rstar+ T_M))))
        ((< H 84.852d3) (let* ((b (b-at-H H))
                               (P (P-at-H-lo-b H b))
                               (T_M (T_M-at-H-lo-b H b)))
                          (/ (* P +M_0+) (* +Rstar+ T_M)))) ;same actually
        (t (/ (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                   (ussadata-v-sumMn-hi *ussadata*)
                           (Z-at-H H))
              +N_A+))))

(defun rho/rho_0-at-H (H)
  "Return density at geopotential altitude, H, in meters divided by the density at altitude zero."
  (/ (mu-at-H H) +mu_0+))

(defun v_m-at-H (H)
  "Return mole volume in m^3 at geopotential altitude, H, in meters."
  (cond ((< H 79.006d3) (let* ((b (b-at-H H))
                               (P (P-at-H-lo-b H b))
                               (Tp (T_M-at-H-lo-b H b))) ;T = T_M in this region
                          (/ (* +Rstar+ Tp) P)))
        ((< H 84.852d3) (let* ((b (b-at-H H))
                               (P (P-at-H-lo-b H b))
                               (Tp (T-at-H-lo-b H b))) ;T != T_M in this region
                          (/ (* +Rstar+ Tp) P)))
        (t (/ +N_A+ (interpolate-uniform (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                                         (ussadata-v-n-hi *ussadata*)
                                 (Z-at-H H))))))

(defun H_P-at-H (H)
  "Return pressure scale height in meters at geopotential altitude, H, in meters."
  (cond ((< H 84.852d3) (let* ((b (b-at-H H))
                               (T_M (T_M-at-H-lo-b H b))
                               (g (g-at-H H)))
                          (/ (* +Rstar+ T_M) (* g +M_0+))))
        (t (let* ((Z (Z-at-H H))
                  (Tp (T-at-Z-hi Z))
                  (n (interpolate-uniform
                      (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                      (ussadata-v-n-hi *ussadata*) Z))
                  (sumMn (interpolate-uniform
                          (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                          (ussadata-v-sumMn-hi *ussadata*) Z))
                  (g (g-at-Z Z)))
             ;; does not quite match the tables in the standard, though its individual components do(?)
             ;; for example at 10000 m, T_M = T = 223.25209264797857d0, g = 9.775868442887434d0,
             ;; and the constants +Rstar+ = 8314.32d0 and +M_0+ = 28.9644d0;
             ;; this gives H_P = 6555.448184469401d0 but the standard gives 6552.5 ??
             (/ (* +Rstar+ Tp n) (* g sumMn))))))

(defun H_rho-at-H (H &key (dZ 100.0d0))
  "Return density scale height in meters at geopotential altitude, H, in meters.  In the high altitude
region, above 86 km, this uses a numerical derivative to get d ln(T_M) / d Z.  DZ is the half-height
over which this derivative is done, centered around Z."
  ;; may as well just do numerical differentiation of ln(rho) then !
  (cond ((< H 84.852d3) (let* ((H_P (H_P-at-H H))
                               (dT/dZ (dT/dZ-at-Z (Z-at-H H)))
                               (Tp (T-at-H-lo H))
                               (dlnTdZ (/ dT/dZ Tp)))
                          (/ H_P (+ 1.0d0 (* H_P dlnTdZ)))))
        (t (let* ((Z (Z-at-H H))
                  (H_P (H_P-at-Z Z))
                  (T_M- (T_M-at-Z (- Z dZ)))
                  (T_M+ (T_M-at-Z (+ Z dZ)))
                  (dlnT_MdZ (/ (- (log T_M+) (log T_M-)) (* 2.0d0 dZ))))
             (/ H_P (+ 1.0d0 (* H_P dlnT_MdZ)))))))

(defun V-at-H (H)
  "Return mean air particle speed in m/s at geopotential altitude, H, in meters."
  (cond ((< H 84.852d3) (let* ((T_M (T_M-at-H-lo H)))
                          (sqrt (/ (* 8.0d0 +Rstar+ T_M) pi +M_0+))))
        (t (let* ((Z (Z-at-H H))
                  (Tp (T-at-Z-hi Z))
                  (n (interpolate-uniform
                      (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                      (ussadata-v-n-hi *ussadata*) Z))
                  (sumMn (interpolate-uniform
                          (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                          (ussadata-v-sumMn-hi *ussadata*) Z)))
             (sqrt (/ (* 8.0d0 +Rstar+ Tp n) pi sumMn))))))

(defun L-at-H (H)
  "Return mean free path in meters at geopotential altitude, H, in meters."
  (cond ((< H 79.006d3) (let* ((b (b-at-H H))
                               (T_M (T_M-at-H-lo-b H b))
                               (P (P-at-H-lo-b H b)))
                          (/ (* +Rstar+ T_M) (* (sqrt 2.0d0) pi +N_A+ +sigma+-sq P))))
        ((< H 84.852d3) (let* ((b (b-at-H H))
                               (T_M (T_M-at-H-lo-b H b))
                               (M/M_0 (M/M_0-at-H-near-86km H))
                               (P (P-at-H-lo-b H b)))
                          (/ (* +Rstar+ M/M_0 T_M)
                             (* (sqrt 2.0d0) pi +N_A+ +sigma+-sq P))))
        (t (let ((n (interpolate-uniform
                     (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                     (ussadata-v-n-hi *ussadata*)
                                 (Z-at-H H))))
             (/ (* (sqrt 2.0d0) pi +sigma+-sq n))))))

(defun nu-at-H (H)
  "Return mean collision frequency in /s at geopotential altitude, H, in meters."
  (cond ((< H 79.006d3) (let* ((b (b-at-H H))
                               (T_M (T_M-at-H-lo-b H b))
                               (P (P-at-H-lo-b H b)))
                          (* 4.0d0 +N_A+ +sigma+-sq
                             (sqrt (/ (* pi P P) (* +Rstar+ +M_0+ T_M))))))
        ((< H 84.852d3) (let* ((b (b-at-H H))
                               (T_M (T_M-at-H-lo-b H b))
                               (M/M_0 (M/M_0-at-H-near-86km H))
                               (M (* +M_0+ M/M_0))
                               (P (P-at-H-lo-b H b)))
                          (* 4.0d0 +N_A+ +sigma+-sq
                             (sqrt (/ (* pi P P) (* +Rstar+ M T_M))))))
        (t (let* ((Z (Z-at-H H))
                  (n (interpolate-uniform
                      (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                      (ussadata-v-n-hi *ussadata*) Z))
                  (sumMn (interpolate-uniform
                          (ussadata-Z-start-hi *ussadata*) (ussadata-dZ-hi *ussadata*)
                          (ussadata-v-sumMn-hi *ussadata*) Z))
                  (P (P-at-Z Z))
                  (Tp (T-at-Z-hi Z)))
             (* 4.0d0 +N_A+ +sigma+-sq
                (sqrt (/ (* pi P P n) (* +Rstar+ Tp sumMn))))))))

(defun C_s-at-H (H)
  "Return speed of sound in m/s at geopotential altitude, H, in meters."
  (sqrt (/ (* +gamma+ +Rstar+ (T_M-at-H H)) +M_0+)))

(defun mu-at-H (H)
  "Return dynamic viscosity in N.s/m^2 at geopotential altitude, H, in meters."
  (let ((Tp (T-at-H H)))
    (/ (* +beta+ (expt Tp 1.5d0)) (+ Tp +S+))))

(defun mu/mu_0-at-H (H)
  "Return dynamic viscosity at geopotential altitude, H, in meters divided by the
dynamic viscosity at altitude zero."
  (/ (mu-at-H H) +mu_0+))

(defun eta-at-H (H)
  "Return kinematic viscosity in m^2/s at geopotential altitude, H, in meters."
  (/ (mu-at-H H) (rho-at-H H)))

(defun eta/eta_0-at-H (H)
  "Return kinematic viscosity at geopotential altitude, H, in meters divided by the
kinematic viscosity at altitude zero."
  (/ (eta-at-H H) +eta_0+))

(defun k_t-at-H (H)
  "Return thermal conductivity in W/m.K at geopotential altitude, H, in meters."
  (let ((Tp (T-at-H H)))
    (/ (* 2.65019d-3 (expt Tp 1.5d0))
       (+ Tp (* 245.4d0 (expt 10.0d0 (/ -12.0d0 Tp)))))))

(defun k_t/k_t0-at-H (H)
  "Return thermal conductivity at geopotential altitude, H, in meters divided by the
thermal conductivity at altitude zero."
  (/ (k_t-at-H H) +k_t0+))


;;;; UNIT CONVERSIONS

;; these are the ones needed to print tables like those in the standard by the functions in
;; test/test-cl-ussa1976.lisp

(defun K->C (K)
  "Return temperature in Celsius corresponding to K in Kelvins."
  (- K 273.15))

(defun Pa->mb (P)
  "Return pressure in millibars corresponding to P in Pascals."
  (* 0.01d0 P))

(defun Pa->torr (P)
  "Return pressure in torr corresponding to P in Pascals."
  (* (/ 760.0d0 101325.0d0) P))
