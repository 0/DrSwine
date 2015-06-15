#lang racket

(require (only-in (planet williams/science/constants)
                  mks-mass-electron
                  mks-molar-gas
                  mks-plancks-constant-hbar
                  mks-unified-atomic-mass
                  num-avogadro))

(provide hbar
         kB
         m-e)


(define hbar (* 1e9 mks-plancks-constant-hbar num-avogadro)) ; kJ ps/mol
(define kB (* 1e-3 mks-molar-gas)) ; kJ/mol K
(define m-e (/ mks-mass-electron mks-unified-atomic-mass)) ; g/mol
