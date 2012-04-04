;; -*- fill-column:102; comment-column:24; -*-

(in-package :test-cl-ussa1976)

(declaim (optimize (debug 3)))
;;(declaim (optimize (speed 3) (safety 0)))

(defun make-table-I-z (&key
                    (path #P"table-I-z.txt")
                    (l-Z (concatenate 'list
                                      (loop for x from -5000.0d0 below 10999.0d0 by 50.0d0
                                         collect x)
                                      (loop for x from 11000.0d0 below 31999.0d0 by 100.0d0
                                           collect x)
                                      (loop for x from 32000.0d0 below 49999.0d0 by 200.0d0
                                           collect x)
                                      (loop for x from 50000.0d0 below 99999.0d0 by 500.0d0
                                           collect x)
                                      (loop for x from 100000.0d0 below 299999.0d0 by 1000.0d0
                                         collect x)
                                      (loop for x from 300000.0d0 below 499999.0d0 by 2000.0d0
                                           collect x)
                                      (loop for x from 500000.0d0 to 1000001.0d0 by 5000.0d0
                                           collect x))))
  "Write a file on PATH with the data in Table I of the standard (pages 50-73)
    Z H T-in-K T-in-C T_M P-in-mb P-in-torr P/P_0 rho rh/rho_0
at the geometric altitudes in meters given by the list L-Z."
  (with-open-file (f path :direction :output :if-exists :supersede :if-does-not-exist :create)
    (format f "~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A~%"
            "Z [m]" "H [m]" "T [K]" "T [C]" "T_M [K]" "P [mb]" "P [torr]" "P/P_0"
            "rho [kg/m^3]" "rho/rho_0" )
    (dolist (Z l-Z)
      (format f
              "~12,0F ~12,0F ~12,3F ~12,3F ~12,3F ~12,4E ~12,4E ~12,4E ~12,4E ~12,4E~%"
              Z (H-at-Z Z) (T-at-Z Z) (K->C (T-at-Z Z)) (T_M-at-Z Z)
              (Pa->mb (P-at-Z Z)) (Pa->torr (P-at-Z Z)) (P/P_0-at-Z Z)
              (rho-at-Z Z) (rho/rho_0-at-Z Z)))))

(defun make-table-I-h (&key
                       (path #P"table-I-h.txt")
                       (l-H (concatenate 'list
                                         (loop for x from -5000.0d0 below 10999.0d0 by 50.0d0
                                            collect x)
                                         (loop for x from 11000.0d0 below 31999.0d0 by 100.0d0
                                            collect x)
                                         (loop for x from 32000.0d0 below 49999.0d0 by 200.0d0
                                            collect x)
                                         (loop for x from 50000.0d0 below 84999.0d0 by 500.0d0
                                            collect x))))
  "Write a file on PATH with the data in Table I of the standard (pages 50-73)
    H Z T-in-K T-in-C T_M P-in-mb P-in-torr P/P_0 rho rh/rho_0
at the geopotential altitudes in meters given by the list L-H."
  (with-open-file (f path :direction :output :if-exists :supersede :if-does-not-exist :create)
    (format f "~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A~%"
            "H [m]" "Z [m]" "T [K]" "T [C]" "T_M [K]" "P [mb]" "P [torr]" "P/P_0"
            "rho [kg/m^3]" "rho/rho_0" )
    (dolist (H l-H)
      (format f
              "~12,0F ~12,0F ~12,3F ~12,3F ~12,3F ~12,4E ~12,4E ~12,4E ~12,4E ~12,4E~%"
              H (Z-at-H H) (T-at-H H) (K->C (T-at-H H)) (T_M-at-H H)
              (Pa->mb (P-at-H H)) (Pa->torr (P-at-H H)) (P/P_0-at-H H)
              (rho-at-H H) (rho/rho_0-at-H H)))))

(defun make-table-II-z (&key
                    (path #P"table-II-z.txt")
                    (l-Z (concatenate 'list
                                      (loop for x from -5000.0d0 below 10999.0d0 by 50.0d0
                                         collect x)
                                      (loop for x from 11000.0d0 below 31999.0d0 by 100.0d0
                                         collect x)
                                      (loop for x from 32000.0d0 below 49999.0d0 by 200.0d0
                                         collect x)
                                      (loop for x from 50000.0d0 below 99999.0d0 by 500.0d0
                                         collect x)
                                      (loop for x from 100000.0d0 below 299999.0d0 by 1000.0d0
                                         collect x)
                                      (loop for x from 300000.0d0 below 499999.0d0 by 2000.0d0
                                         collect x)
                                      (loop for x from 500000.0d0 to 1000001.0d0 by 5000.0d0
                                         collect x))))
  "Write a file on PATH with the data in Table II of the standard (starting page 74-97)
    Z H g H_p n V nu L M
at the geometric altitudes in meters given by the list L-Z."
  (with-open-file (f path :direction :output :if-exists :supersede :if-does-not-exist :create)
    (format f "~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A~%"
            "Z" "H" "g" "H_p" "n" "V" "nu" "L" "M")
    (dolist (Z l-Z)
      (format f
              "~12,0F ~12,0F ~12,4F ~12,1F ~12,4E ~12,2F ~12,4E ~12,4E ~12,3F~%"
              Z (H-at-Z Z) (g-at-Z Z) (H_p-at-Z Z) (n-at-Z Z) (V-at-Z Z) (nu-at-Z Z)
              (L-at-Z Z) (M-at-Z Z)))))

(defun make-table-II-h (&key
                        (path #P"table-II-h.txt")
                        (l-H (concatenate 'list
                                          (loop for x from -5000.0d0 below 10999.0d0 by 50.0d0
                                             collect x)
                                          (loop for x from 11000.0d0 below 31999.0d0 by 100.0d0
                                             collect x)
                                          (loop for x from 32000.0d0 below 49999.0d0 by 200.0d0
                                             collect x)
                                          (loop for x from 50000.0d0 below 84999.0d0 by 500.0d0
                                             collect x))))
  "Write a file on PATH with the data in Table II of the standard (starting page 74-97)
    Z H g H_p n V nu L M
at the geopotential altitudes in meters given by the list L-H."
  (with-open-file (f path :direction :output :if-exists :supersede :if-does-not-exist :create)
    (format f "~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A~%"
            "Z" "H" "g" "H_p" "n" "V" "nu" "L" "M")
    (dolist (H l-H)
      (format f
              "~12,0F ~12,0F ~12,4F ~12,1F ~12,4E ~12,2F ~12,4E ~12,4E ~12,3F~%"
              H (Z-at-H H) (g-at-H H) (H_p-at-H H) (n-at-H H) (V-at-H H) (nu-at-H H)
              (L-at-H H) (M-at-H H)))))

(defun make-table-III-z (&key
                    (path #P"table-III-z.txt")
                    (l-Z (concatenate 'list
                                      (loop for x from -5000.0d0 below 10999.0d0 by 50.0d0
                                         collect x)
                                      (loop for x from 11000.0d0 below 31999.0d0 by 100.0d0
                                         collect x)
                                      (loop for x from 32000.0d0 below 49999.0d0 by 200.0d0
                                         collect x)
                                      (loop for x from 50000.0d0 below 85599.0d0 by 500.0d0
                                         collect x))))
  "Write a file on PATH with the data in Table III of the standard (pages 98-115)
    Z H C_s mu mu/mu_0 eta eta/eta_0 k_t k_t/k_t0
at the geometric altitudes in meters given by the list L-Z."
  (with-open-file (f path :direction :output :if-exists :supersede :if-does-not-exist :create)
    (format f "~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A~%"
            "Z" "H" "C_s" "mu" "mu/mu_0" "eta" "eta/eta_0" "k_t" "k_t/k_t0")
    (dolist (Z l-Z)
      (format f
              "~12,0F ~12,0F ~12,2F ~12,4E ~12,4E ~12,4E ~12,4E ~12,4E ~12,4E~%"
              Z (H-at-Z Z) (C_s-at-Z Z) (mu-at-Z Z) (mu/mu_0-at-Z Z)
              (eta-at-Z Z) (eta/eta_0-at-Z Z) (k_t-at-Z Z) (k_t/k_t0-at-Z Z)))))

(defun make-table-III-h (&key
                         (path #P"table-III-h.txt")
                         (l-H (concatenate 'list
                                           (loop for x from -5000.0d0 below 10999.0d0 by 50.0d0
                                              collect x)
                                           (loop for x from 11000.0d0 below 31999.0d0 by 100.0d0
                                              collect x)
                                           (loop for x from 32000.0d0 below 49999.0d0 by 200.0d0
                                              collect x)
                                           (loop for x from 50000.0d0 below 84599.0d0 by 500.0d0
                                              collect x))))
  "Write a file on PATH with the data in Table III of the standard (pages 98-115)
    H Z C_s mu mu/mu_0 eta eta/eta_0 k_t k_t/k_t0
at the geopotential altitudes in meters given by the list L-H."
  (with-open-file (f path :direction :output :if-exists :supersede :if-does-not-exist :create)
    (format f "~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A~%"
            "H" "Z" "C_s" "mu" "mu/mu_0" "eta" "eta/eta_0" "k_t" "k_t/k_t0")
    (dolist (H l-H)
      (format f
              "~12,0F ~12,0F ~12,2F ~12,4E ~12,4E ~12,4E ~12,4E ~12,4E ~12,4E~%"
              H (Z-at-H H) (C_s-at-H H) (mu-at-H H) (mu/mu_0-at-H H)
              (eta-at-H H) (eta/eta_0-at-H H) (k_t-at-H H) (k_t/k_t0-at-H H)))))

(defun make-table-VIII-z (&key
                    (path #P"table-VIII-z.txt")
                    (l-Z (concatenate 'list
                                      (loop for x from 86000.0d0 below 99999.0d0 by 500.0d0
                                         collect x)
                                      (loop for x from 100000.0d0 below 299999.0d0 by 1000.0d0
                                         collect x)
                                      (loop for x from 300000.0d0 below 499999.0d0 by 2000.0d0
                                         collect x)
                                      (loop for x from 500000.0d0 to 1000001.0d0 by 5000.0d0
                                         collect x))))
  "Write a file on PATH with the data in Table VIII of the standard (pages 210-215)
    number density of individual species, N2 O O2 Ar He H
    and total number density (in addition to the individual species data which the standard has)
at the geometric altitudes in meters given by the list L-Z."
  (with-open-file (f path :direction :output :if-exists :supersede :if-does-not-exist :create)
    (format f "~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A ~12@A~%"
            "Z" "H" "n_N2" "n_O" "n_O2" "n_Ar" "n_He" "n_H" "n")
    (dolist (Z l-Z)
      (let* ((v-n (v-n_X-at-Z Z))
             (n_N2 (aref v-n 0))
             (n_O (aref v-n 1))
             (n_O2 (aref v-n 2))
             (n_Ar (aref v-n 3))
             (n_He (aref v-n 4))
             (n_H (aref v-n 5))
             (n (aref v-n 6)))
        (format f
                "~12,0F ~12,0F ~12,3E ~12,3E ~12,3E ~12,3E ~12,3E ~12,3E   ~12,4E~%"
                Z (H-at-Z Z) n_N2 n_O n_O2 n_Ar n_He n_H n)))))
