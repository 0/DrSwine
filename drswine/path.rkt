#lang racket

(require (for-syntax
           racket/syntax))

(require (only-in racket/class
                  class*))

(require (only-in math/distributions
                  normal-dist
                  sample))

(require (only-in "constants.rkt"
                  hbar
                  kB)
         (only-in "fft.rkt"
                  dct
                  idct))

(provide path%)


(define (random-normal [var 1.0] #:num [num #f])
  (let* ([d (normal-dist 0.0 (sqrt var))]
         [n (or num 1)]
         [s (sample d n)])
    (if num (list->vector s) (car s))))

(define (random-momenta mass temperature #:num [num #f])
  (random-normal (* mass kB temperature) #:num num))

; Automatically convert momenta and positions to and from normal modes as
; necessary.
(define-syntax (auto-nm stx)
  (syntax-case stx ()
    [(auto-nm to x)
     (with-syntax ([xs (datum->syntax stx (format-symbol "_~as" #'x))]
                   [nmxs (datum->syntax stx (format-symbol "_nm~as" #'x))]
                   [in-nmx (datum->syntax stx (format-symbol "_in-nm~a" #'x))])
       (let* ([to_ (syntax->datum #'to)]
              [this (if to_ #'nmxs #'xs)])
         #`(if #,(if to_ #'in-nmx #'(not in-nmx))
             #,this
             (begin
               (set! #,this (#,(if to_ 'dct 'idct) #,(if to_ #'xs #'nmxs)))
               (set! in-nmx to)
               #,this))))]))

(define path%
  (class* object% (printable<%>)
          (init num-beads
                bead-mass
                tau
                [centroid-friction 0.0]
                [name #f])
          (super-new)

          (define _bead-mass bead-mass)
          (define/public (get-bead-mass) _bead-mass)

          (define _tau tau)
          (define/public (get-tau) _tau)

          (define _num-beads num-beads)
          (define/public (get-num-beads) _num-beads)

          (define _centroid-friction centroid-friction)
          (define/public (get-centroid-friction) _centroid-friction)

          (define _name name)
          (define/public (get-name) _name)

          ; Start in Cartesian coordinates.
          (define _in-nmp #f)
          (define _in-nmq #f)

          (define _ps (make-vector _num-beads)) ; g nm/ps mol
          (define/public (get-ps) (auto-nm #f p))

          (define _qs (make-vector _num-beads)) ; nm
          (define/public (get-qs) (auto-nm #f q))

          (define _nmps (void))
          (define/public (get-nmps) (auto-nm #t p))

          (define _nmqs (void))
          (define/public (get-nmqs) (auto-nm #t q))

          (define _omega-ks (build-vector _num-beads
                                          (lambda (k) (/ (sin (/ (* 0.5 pi k)
                                                                 _num-beads))
                                                         (* 0.5 hbar _tau)))))
          (define _gammas (vector-map (lambda (omega-k)
                                        (if (= 0.0 omega-k)
                                          _centroid-friction
                                          (* 2.0 omega-k)))
                                      _omega-ks))

          (define/public (initialize-momenta-to-temperature temperature)
                         (vector-map! (lambda (_ p) p)
                                      (get-ps)
                                      (random-momenta _bead-mass
                                                      temperature
                                                      #:num _num-beads))
                         (void))

          (define/public (kinetic-energy) ; kJ/mol
                         (/ (for/sum ([p (get-ps)]) (* p p))
                            (* 2.0 _bead-mass)))

          (define/public (spring-energy) ; kJ/mol
                         (for/sum ([nmq (get-nmqs)]
                                   [omega-k _omega-ks])
                                  (* 0.5 _bead-mass (sqr omega-k) (sqr nmq))))

          (define/public (propagate dt)
                         (define cs (vector-map (lambda (omega-k)
                                                  (cos (* omega-k dt)))
                                                _omega-ks))
                         (define ss (vector-map (lambda (omega-k)
                                                  (sin (* omega-k dt)))
                                                _omega-ks))
                         (define sds (vector-map (lambda (s omega-k)
                                                   (if (= 0.0 omega-k)
                                                     dt
                                                     (/ s omega-k)))
                                                 ss _omega-ks))
                         (define temp-nmps (vector-copy (get-nmps)))
                         (vector-map! (lambda (p q omega-k c s)
                                        (- (* c p)
                                           (* _bead-mass omega-k s q)))
                                      (get-nmps) (get-nmqs) _omega-ks cs ss)
                         (vector-map! (lambda (q p omega-k c sd)
                                        (+ (* c q)
                                           (/ (* sd p) _bead-mass)))
                                      (get-nmqs) temp-nmps _omega-ks cs sds)
                         (void))

          (define/public (thermostat temperature dt)
                         (vector-map! (lambda (nmp gamma r)
                                        (let* ([c1 (exp (* -0.5 gamma dt))]
                                               [c2 (sqrt (* _bead-mass kB temperature (- 1.0 (sqr c1))))])
                                          (+ (* c1 nmp)
                                             (* c2 r))))
                                      (get-nmps)
                                      _gammas
                                      (random-normal #:num _num-beads))
                         (void))

          (define/public (custom-print p depth) (display _name p))
          (define/public (custom-write p)       (write _name p))
          (define/public (custom-display p)     (display _name p))))
