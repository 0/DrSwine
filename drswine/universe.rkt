#lang racket

(require (only-in racket/class
                  class*))

(require (only-in srfi/26
                  cut))

(require (only-in "constants.rkt"
                  kB))

(provide add-path
         default-universe
         initialize-all-momenta
         universe%)


(define universe%
  (class object%
         (init temperature)
         (super-new)

         (define _temperature temperature)
         (define/public (get-temperature) _temperature)

         (define _paths '())
         (define/public (get-paths) _paths)
         (define/public (add-path p)
                        (set! _paths (cons p _paths)))
         (define/public (find-path name)
                        (car (memf (lambda (p) (string=? (send p get-name) name))
                                   _paths)))
         (define/public (remove-path p)
                        (set! _paths (remove p _paths)))

         (define/public (initialize-all-momenta)
                        (for-each (cut send <> initialize-momenta-to-temperature _temperature)
                                  _paths))

         (define/public (num-beads)
                        (for/sum ([p _paths]) (send p get-num-beads)))

         (define/public (kinetic-energy) ; kJ/mol
                        (for/sum ([p _paths]) (send p kinetic-energy)))

         (define/public (spring-energy) ; kJ/mol
                        (for/sum ([p _paths]) (send p spring-energy)))

         (define/public (effective-temperature) ; K
                        (/ (* 2.0 (kinetic-energy))
                           (* kB (num-beads))))))


(define default-universe
  (make-parameter
    #f
    (lambda (u)
      (if (is-a? u universe%)
        u
        (error "not a universe%" u)))))


(define (add-path path [universe (default-universe)])
  (unless universe
    (error "no universe specified"))
  (send universe add-path path))

(define (initialize-all-momenta [universe (default-universe)])
  (unless universe
    (error "no universe specified"))
  (send universe initialize-all-momenta))
