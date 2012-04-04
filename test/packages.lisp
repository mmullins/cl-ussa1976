;; -*- fill-column:102; comment-column:24; -*-

(defpackage :test-cl-ussa1976
  (:use :common-lisp :cl-ussa1976)
  (:documentation
   "Test Common Lisp implementation of US Standard Atmosphere, 1976.

Function MAKE-TABLE-I-Z writes onto a file this implementation's values for the data in Table I
of the standard, containing
    Z H T-in-K T-in-C T_M P-in-mb P-in-torr P/P_0 rho rh/rho_0
at specified geometric altitudes, Z.

Function MAKE-TABLE-II-Z writes onto a file this implementation's values for the data in Table II
of the standard, containing
    Z H g H_p n V nu L M
at specified geometric altitudes, Z.

Function MAKE-TABLE-III-Z writes onto a file this implementation's values for the data in Table III
of the standard, containing
    Z H C_s mu mu/mu_0 eta eta/eta_0 k_t k_t/k_t0
at specified geometric altitudes, Z.

Function MAKE-TABLE-VIII-Z writes onto a file this implementation's values for the data in Table VIII
of the standard, containing
    number density of individual species, N2 O O2 Ar He H
    and total number density (in addition to the individual species data which the standard has)
at specified geometric altitudes, Z.

Functions MAKE-TABLE-I-H, MAKE-TABLE-II-H, MAKE-TABLE-III-H are analogous to the ...-Z functions
except that they take geopotential altitudes, H, instead of geometric altitudes, Z.  The standard
gives data in that format for these tables for H values below 84500 m.

Note that in a number of cases in the tables in the standard, the headers of columns Z and H are
switched.  Also in Table II, the exponent of thermal conductivity, k, should be -2 not -5.  See the
corrigenda at the beginning of the standard for details.")
  (:export
   #:make-table-I-z #:make-table-I-h
   #:make-table-II-z #:make-table-II-h
   #:make-table-III-z #:make-table-III-h
   #:make-table-VIII-z))
