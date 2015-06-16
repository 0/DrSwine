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

         (define _next-bead-index 0)
         (define _bead-index-mapping (make-hash))

         (define _paths '())
         (define/public (get-paths) _paths)
         (define/public (add-path p)
                        (send p set-universe this)
                        (let* ([num-beads (send p get-num-beads)]
                               [existing-indices (send p get-bead-indices)]
                               [indices
                                 (if existing-indices
                                   existing-indices
                                   (build-vector num-beads
                                                 (cut + <> _next-bead-index)))])
                          (unless existing-indices
                            (send p set-bead-indices indices)
                            (set! _next-bead-index (+ num-beads _next-bead-index)))
                          (for ([i (in-range num-beads)]
                                [idx indices])
                            (hash-set! _bead-index-mapping idx (cons p i))))
                        (set! _paths (cons p _paths)))
         (define/public (find-path name)
                        (car (memf (lambda (p) (string=? (send p get-name) name))
                                   _paths)))
         (define/public (remove-path p)
                        (set! _paths (remove p _paths))
                        (for ([idx (send p get-bead-indices)])
                          (hash-remove! _bead-index-mapping idx))
                        (send p clear-universe))

         (define/public (initialize-all-momenta)
                        (for-each (cut send <> initialize-momenta-to-temperature _temperature)
                                  _paths))

         (define/public (num-beads)
                        (for/sum ([p _paths]) (send p get-num-beads)))

         (define/public (kinetic-energy) ; kJ/mol
                        (for/sum ([p _paths]) (send p kinetic-energy)))

         (define/public (spring-energy) ; kJ/mol
                        (for/sum ([p _paths]) (send p spring-energy)))

         (define _potential-energy #f)
         (define/public (potential-energy) ; kJ/mol
                        _potential-energy)
         (define/public (set-potential-energy energy)
                        (set! _potential-energy energy))

         (define/public (effective-temperature) ; K
                        (/ (* 2.0 (kinetic-energy))
                           (* kB (num-beads))))

         (define/public (get-position idx) ; nm
                        (match-let ([(cons p i) (hash-ref _bead-index-mapping idx)])
                          (vector-ref (send p get-qs) i)))

         (define/public (apply-force idx force)
                        (match-let ([(cons p i) (hash-ref _bead-index-mapping idx)])
                          (let ([ps (send p get-ps)])
                            (vector-set! ps i
                                         (+ (vector-ref ps i)
                                            force)))))))


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
