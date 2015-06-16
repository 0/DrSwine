#lang racket

#|
 | Sample from the distribution of middle bead positions for a harmonic
 | oscillator.
 |
 | Suggested parameters:
 |
 |   --mass 1.0 --omega 1.0 --beta 1.0 --num-links 8 --dt 0.1 --num-steps 500000 --trial-deform 1.1
 |#

(require (only-in "../drswine/constants.rkt"
                  hbar
                  kB
                  m-e)
         (only-in "../drswine/force-field.rkt"
                  apply-ff-to-paths
                  ff-harmonic-trap)
         (only-in "../drswine/integrator.rkt"
                  integrator%)
         (only-in "../drswine/path.rkt"
                  path%)
         (only-in "../drswine/universe.rkt"
                  add-path
                  default-universe
                  initialize-all-momenta
                  universe%))


;;; Parse arguments.
(define cmd-req-params '())

(define-syntax-rule
  (define-cmd-req-param label name)
  (begin
    (define name (make-parameter #f))
    (set! cmd-req-params (cons (cons label name) cmd-req-params))))

(define-cmd-req-param "mass" par-mass-e)
(define-cmd-req-param "omega" par-omega-K)
(define-cmd-req-param "beta" par-beta-K)
(define-cmd-req-param "num-links" par-num-links)
(define-cmd-req-param "dt" par-dt)
(define-cmd-req-param "num-steps" par-num-steps)
(define par-centroid-friction (make-parameter 0.0))
(define par-trial-deform (make-parameter 1.0))
(define par-output-data (make-parameter #f))

(command-line
  #:once-each
  [("--mass") M "Mass (electron masses)"
              (par-mass-e (string->number M))]
  [("--omega") W "Angular frequency (K)"
              (par-omega-K (string->number W))]
  [("--beta") B "Beta (1/K)"
              (par-beta-K (string->number B))]
  [("--num-links") N "Number of links"
                   (par-num-links (string->number N))]
  [("--dt") T "Time step (ps)"
            (par-dt (string->number T))]
  [("--num-steps") N "Number of steps"
                   (par-num-steps (string->number N))]
  [("--centroid-friction") G "Centroid friction (1/ps; default: 0.0)"
                           (par-centroid-friction (string->number G))]
  [("--trial-deform") D "Trial function deformation (default: 1.0)"
                           (par-trial-deform (string->number D))]
  [("--output-data") F "Output path (default: stdout)"
                     (par-output-data F)])

(for ([x cmd-req-params])
  (unless ((cdr x))
    (printf "~a not given~n" (car x))
    (exit)))

(define mass (* (par-mass-e) m-e)) ; g/mol
(define omega (/ (* (par-omega-K) kB) hbar)) ; 1/ps
(define beta-K (par-beta-K)) ; 1/K
(define num-links (par-num-links)) ; 1
(define dt (par-dt)) ; ps
(define num-steps (par-num-steps)) ; 1
(define centroid-friction (par-centroid-friction)) ; 1/ps
(define trial-deform (par-trial-deform)) ; 1

(unless (even? num-links)
  (error "number of links must be even"))

(define bead-mass (/ mass num-links)) ; g/mol
(define temperature (/ 1.0 beta-K)) ; K
(define beta (/ beta-K kB)) ; mol/kJ
(define tau (/ beta num-links)) ; mol/kJ
(define num-beads (add1 num-links)) ; 1

; Output exact distribution parameters.
(parameterize ([current-output-port (current-error-port)])
  (let ([var (/ hbar (* 2.0 mass omega))])
    (printf "variance: ~a nm^2~n" var)
    (printf "normalization: ~a 1/nm~n" (sqrt (/ (* 2.0 pi var))))))


;;; Make objects.
(define u (new universe%
               [temperature temperature]))
(default-universe u)

(define p1 (new path%
                [num-beads num-beads]
                [bead-mass bead-mass]
                [tau tau]
                [centroid-friction centroid-friction]
                [name "p1"]))
(add-path p1)


(define force-constant (* mass (sqr omega))) ; kJ/mol nm^2
(define force-constant-trial (/ (* trial-deform mass omega)
                                hbar)) ; 1/nm^2
(define scaling-trial (/ 1.0 beta)) ; kJ/mol

(define force-fields
  (list
    (apply-ff-to-paths
      (ff-harmonic-trap force-constant 0.0)
      p1
      #:scale-all 'links
      #:scale-ends 'pigs
      #:name "harmonic-trap")
    (apply-ff-to-paths
      (ff-harmonic-trap force-constant-trial 0.0)
      p1
      #:scale-ends scaling-trial
      #:only-ends #t
      #:name "trial")))


(define itg (new integrator%
                 [universe u]
                 [dt dt]
                 [force-fields force-fields]))


;;; Equilibrate.

(let ([runs 10]
      [steps 1000])
  (for ([_ (in-range runs)])
    (initialize-all-momenta)
    (for ([_ (in-range steps)])
      (send itg step))))

(initialize-all-momenta)


;;; Production run.

(define (step-callback _)
  (let ([qs (send p1 get-qs)]
        [M (quotient num-beads 2)])
    (displayln (vector-ref qs M))))

(define (run)
  (for ([_ (in-range num-steps)])
    (send itg step step-callback)))

(if (par-output-data)
    (with-output-to-file
      (par-output-data) #:mode 'text #:exists 'truncate/replace
      run)
    (run))
