;; -*- fill-column:102; comment-column:24; -*-

(defpackage :cl-ussa1976
  (:nicknames :ussa)
  (:documentation
   "US Standard Atmosphere, 1976, per 'US Standard Atmosphere, 1976', NOAA-S/T 76-1562.
Units are [m] [kg] [s] [K].")
  (:use :common-lisp)
  (:export

   ;; USSADATA STRUCTURE
   #:make-ussadata-with-dZ-hi
   #:with-ussadata

   ;; PROPERTIES AT GEOMETRIC ALTITUDE Z
   #:H-at-Z
   #:g-at-Z
   #:M/M_0-at-Z
   #:M-at-Z
   #:T-at-Z
   #:T_M-at-Z
   #:n-at-Z
   #:n_N2-at-Z
   #:n_O-at-Z
   #:n_O2-at-Z
   #:n_Ar-at-Z
   #:n_He-at-Z
   #:n_H-at-Z
   #:v-n_X-at-Z
   #:P-at-Z
   #:P/P_0-at-Z
   #:rho-at-Z
   #:rho/rho_0-at-Z
   #:v_m-at-Z
   #:H_P-at-Z
   #:H_rho-at-Z
   #:V-at-Z
   #:L-at-Z
   #:nu-at-Z
   #:C_s-at-Z
   #:mu-at-Z
   #:mu/mu_0-at-Z
   #:eta-at-Z
   #:eta/eta_0-at-Z
   #:k_t-at-Z
   #:k_t/k_t0-at-Z

   ;; PROPERTIES AT GEOPOTENTIAL ALTITUDE H
   #:Z-at-H
   #:g-at-H
   #:M/M_0-at-H
   #:M-at-H
   #:T-at-H
   #:T_M-at-H
   #:n-at-H
   #:n_N2-at-H
   #:n_O-at-H
   #:n_O2-at-H
   #:n_Ar-at-H
   #:n_He-at-H
   #:n_H-at-H
   #:v-n_X-at-H
   #:P-at-H
   #:P/P_0-at-H
   #:rho-at-H
   #:rho/rho_0-at-H
   #:v_m-at-H
   #:H_P-at-H
   #:H_rho-at-H
   #:V-at-H
   #:L-at-H
   #:nu-at-H
   #:C_s-at-H
   #:mu-at-H
   #:mu/mu_0-at-H
   #:eta-at-H
   #:eta/eta_0-at-H
   #:k_t-at-H
   #:k_t/k_t0-at-H

   ;;;; UNIT CONVERSIONS
   #:K->C
   #:Pa->mb
   #:Pa->torr))
