isAsmJs= false;

/*
function Test( stdlib, foreign, heap ) {
    "use asm";

    function run() {
        var f= 0.0;
        var e= 0;
        f= e;
    }
    
    return {};
}
*/
  
function Raytracer( stdlib, foreign, heap ) {
//    "use asm";

    var sqrt = stdlib.Math.sqrt;
    var imul = stdlib.Math.imul;
    var pow = stdlib.Math.pow;
    
    // var heap8 = new stdlib.Uint8Array(heap);
    var heap8= heap;

    var w= 0;
    var inx_spheres= 0;
    var ambient_light= 0.0;
    var inx_lights= 0;
    var inx_output= 0;

    var t= 0.0;
    var c= 0;
    
    function max( a, b ) {
        a= +a;
        b= +b;
        return +(a > b ? a : b);
    }

    // Dot product.
    function dot(Ax, Ay, Az, Bx, By, Bz) {
        Ax= +Ax; Ay= +Ay; Az= +Az;
        Bx= +Bx; By= +By; Bz= +Bz;
        return +(Ax * Bx + Ay * By + Az * Bz);
    }

    // -----------------------------------------------------------------------------

    // Shorten some names.
    // var math = Math;
    // var sqrt = math.sqrt;
    // var max = math.max;

    // Closure doesn't rename vars unless they're declared with "var", which takes
    // space. So most vars are 1-letter and global:
    //
    // C: sphere center
    // L: light vector
    // N: surface normal at intersection
    // X: intersection point
    // a: quadratic equation constant
    // b: quadratic equation constant
    // c: color channel
    // d: quadratic equation discriminant
    // e: loop variable
    // f: candidate parameter t
    // h: half-width of the canvas
    // i: illumination
    // J: (ray origin) - (sphere center) 
    // k: <N, L> 
    // l: light index in loop
    // n: <N, N>
    // q: sphere index in loop
    // r: sphere radius
    // s: closest intersection sphere index
    // t: closest intersection t
    // u: intensity of lights[l] 
    // v: closest sphere found in loop
    //
    // The exceptions are vars that need to be initialized here (we still pay the
    // "a=", so we pay a single "var" above, and use nice names) and some vars in 
    // trace_ray, which is recursive, so some of it vars can't be global.

    // Helper: A_minus_Bk(A, B, k)  =  A - B*k. Since it's used more with k < 0,
    // using - here saves a couple of bytes later.
    // function A_minus_Bk (A, B, k) {
    //     return [A[0] - B[0]*k, A[1] - B[1]*k, A[2] - B[2]*k];
    // }


    // Find nearest intersection of the ray from B in direction D with any sphere.
    // "Interesting" parameter values must be in the range [t_min, t_max].
    // Returns the index within spheres of the center of the hit sphere, 0 if none.
    // The parameter value for the intersection is in the global variable t.
    function closest_intersection(Bx, By, Bz, Dx, Dy, Dz, t_min, t_max) {
        Bx= +Bx; By= +By; Bz= +Bz;
        Dx= +Dx; Dy= +Dy; Dz= +Dz;
        t_min= +t_min;
        t_max= +t_max;

        var v= 0;
        var q= 0;
        var r= 0;
        var Jx= 0.0, Jy= 0.0, Jz= 0.0;
        var a= 0.0;
        var b= 0.0;
        var d= 0.0;
        var e= 0;
        var f= 0.0;

        t = +(w|0);  // Min distance found.

        // For each sphere.
        // Get the radius and test for end of array at the same time; 
        // spheres[n] == undefined ends the loop.
        // q points to the 2nd element of the sphere because of q++; +6 skips to next
        // sphere.

        for ( ; r = heap8[inx_spheres + (q|0)]|0; q= (q|0) + 9) {  
            // Compute quadratic equation coefficients K1, K2, K3

            Jx = +(Bx - +(heap8[inx_spheres + q + (1|0)] >> 0));  // origin - center
            Jy = +(By - +(heap8[inx_spheres + q + (2|0)] >> 0));  // origin - center
            Jz = +(Bz - +(heap8[inx_spheres + q + (3|0)] >> 0));  // origin - center
            
            a =  2.0 * +dot(Dx, Dy, Dz, Dx, Dy, Dz);  // 2*K1
            b = -2.0 * +dot(Jx, Jy, Jz, Dx, Dy, Dz);  // -K2
// console.log(r, a, b, Jx, Jy, Jz, Dx, Dy, Dz);
// console.log(By, heap8[inx_spheres + q + (2|0)], (heap8[inx_spheres + q + (2|0)] >>> 0));

            // Compute sqrt(Discriminant) = sqrt(K2*K2 - 4*K1*K3), go ahead if there are
            // solutions.
            d = +sqrt(b * b - 2.0 * a * (+dot(Jx, Jy, Jz, Jx, Jy, Jz) - +((r|0) * (r|0))))
            if ( d != 0.0 ) {
                // Compute the two solutions.
                for (e = 2; e; e = (e|0) - 1, d = -d) {
                    f = (b - d) / a;  // f = (-K2 - d) / 2*K1
// console.log("f", b, b * b, 2.0 * a * +dot(Jx, Jy, Jz, Jx, Jy, Jz), (r|0 * r), "e", e);
                    if ( t_min < f ) if ( f < t_max ) if ( f < t ) { 
                        v = (q|0) + 1;
                        t = f;
                    }
                }
            }
        }

        // Return index of closest sphere in range; t is global
        return v|0;
    }


    // Trace the ray from B with direction D considering hits in [t_min, t_max].
    // If depth > 0, trace recursive reflection rays.
    // Returns the value of the current color channel as "seen" through the ray.
    function trace_ray(Bx, By, Bz, Dx, Dy, Dz, t_min, t_max, depth) {

        var s, Nx, Ny, Nz, Xx, Xy, Xz, n, i, l, u, k, Lx, Ly, Lz, Mx, My, Mz;
// console.log(depth, Dx, Dy, Dz);

        // Find nearest hit; if no hit, return black.
        if (!(s = closest_intersection(Bx, By, Bz, Dx, Dy, Dz, t_min, t_max))) {
            return 0;
        }

        // Compute "normal" at intersection: N = X - spheres[s]
        // intersection: X = B + D*t = B - D*(-t)
        // X = A_minus_Bk(Bx, By, Bz, Dx, Dy, Dz, -t);
        Xx = Bx + Dx * t;
        Xy = By + Dy * t;
        Xz = Bz + Dz * t;
        
        Nx = Xx - heap8[inx_spheres + s];
        Ny = Xy - heap8[inx_spheres + s + 1];
        Nz = Xz - heap8[inx_spheres + s + 2];

        // Instead of normalizing N, we divide by its length when appropriate. Most of
        // the time N appears twice, so we precompute its squared length.
        n = dot(Nx, Ny, Nz, Nx, Ny, Nz);

        // Start with ambient light only
        i = ambient_light;

        // For each light
        for (l = 0; u = heap8[inx_lights + l]; l = l + 4) { // Get intensity and check for end of array

            // Compute vector from intersection to light (L = lights[l++] - X) and
            // k = <N,L> (reused below)
            // L = A_minus_Bk(lights[l++], Xx, Xy, Xz, 1);
            Lx = heap8[inx_lights + l + 1] - Xx;
            Ly = heap8[inx_lights + l + 2] - Xy;
            Lz = heap8[inx_lights + l + 3] - Xz;

            k = dot(Nx, Ny, Nz, Lx, Ly, Lz);

            // Add to lighting

            // If the point isn't in shadow
            // [t_min, t_max]  =  [epsilon,  1] - epsilon avoids self-shadow, 1 
            // doesn't look farther than the light itself.

            // Diffuse lighting, only if it's facing the point 
            // <N,L> / (|N|*|L|) = cos(alpha)
            // Also, |N|*|L| = sqrt(<N,N>)*sqrt(<L,L>) = sqrt(<N,N>*<L,L>)

            // Specular highlights
            //
            // specular = (<R,V>   / (|R|*|V|))   ^ exponent
            //          = (<-R,-V> / (|-R|*|-V|)) ^ exponent
            //          = (<-R,D>  / (|-R|*|D|))  ^ exponent
            //
            // R = 2*N*<N,L> - L
            // M = -R = -2*N*<N,L> + L = L + N*(-2*<N,L>)
            //
            // If the resultant intensity is negative, treat it as 0 (ignore it).

            var cc= closest_intersection(Xx, Xy, Xz, Lx, Ly, Lz, 1 / w, 1);

            // M = A_minus_Bk(Lx, Ly, Lz, Nx, Ny, Nz, 2*k/n), D) 
            var mf= 2 * k / n;
            Mx = Lx - Nx * mf;
            My = Ly - Ny * mf;
            Mz = Lz - Nz * mf;
            
            i += u * 
                !cc * (
                    max(0, k / sqrt(dot(Lx, Ly, Lz, Lx, Ly, Lz) * n))
                    + max(0, pow(dot(Mx, My, Mz, Dx, Dy, Dz) 
                        / sqrt(dot(Mx, My, Mz, Mx, My, Mz) * dot(Dx, Dy, Dz, Dx, Dy, Dz)), heap8[inx_spheres + s + 6]))
                );
        }


        // Compute the color channel multiplied by the light intensity. 2.8 maps
        // the color range from [0, 9] to [0, 255] and the intensity from [0, 10]
        // to [0, 1],  because 2.8 ~ (255/9)/10
        // 
        // spheres[s] = sphere center, so spheres[s+c] = color channel
        // (c = [1..3] because ++c below)
        var local_color = heap8[inx_spheres + s + c + 2] * i * 2.8;

        // If the recursion limit hasn't been hit yet, trace reflection rays.
        // N = normal (non-normalized - two divs by |N| = div by <N,N>
        // D = -view
        // R = 2*N*<N,V>/<N,N> - V = 2*N*<N,-D>/<N,N> + D = D - N*(2*<N,D>/<N,N>)
        var ref = heap8[inx_spheres + s + 7] / 9;
        
        if ( depth-- ) {

            var Rx, Ry, Rz;
            var rf = 2 * dot(Nx, Ny, Nz, Dx, Dy, Dz) / n;
            Rx= Dx - Nx * rf;
            Ry= Dy - Ny * rf;
            Rz= Dz - Nz * rf;

            return trace_ray(Xx, Xy, Xz, Rx, Ry, Rz, 1/w, w, depth) * ref
                + local_color * (1 - ref);
        }
                
        return local_color;
    }
    
    function run( heapIndex ) {
        heapIndex= heapIndex|0;
    
        // Tiny Raytracer (C) Gabriel Gambetta 2013
        // ----------------------------------------
        //
        //  Configuration and scene
        //
        // Size of the canvas. w is also reused as a "big constant" / "+infinity"
        
        w = heap8[heap8[heapIndex]|0]|0;

        // Sphere: radius, [cx,  cy,  cz], R,  G,  B, specular exponent, reflectiveness 
        // R, G, B in [0, 9], reflectiveness in [0..9].
        inx_spheres= heap8[(heapIndex|0) + 1]|0;

        // Ambient light.
        ambient_light = heap8[heap8[(heapIndex|0) + 2]|0]|0;

        // Point lights: intensity, [x,  y,  z]
        // Intensities should add to 10, including ambient.
        inx_lights = heap8[(heapIndex|0) + 3]|0;

        inx_output= heap8[(heapIndex|0) + 4]|0;

        var y, h, x;

        var out_idx = inx_output + 1;

        // For each y; also compute h=w/2 without paying an extra ";"
        for (y = h=w/2; y-- > -h;) {

            // For each x
            for (x = -h; x++ < h;) {

                // One pass per color channel (!). This way we don't have to deal with
                // "colors".
                for (c = 0; ++c < 4;) {
                    // Camera is at (0, 1, 0)
                    //
                    // Ray direction is (x*vw/cw, y*vh/ch, 1) where vw = viewport width, 
                    // cw = canvas width (vh and ch are the same for height). vw is fixed
                    // at 1 so (x/w, y/w, 1)
                    //
                    // [t_min, t_max] = [1, w], 1 starts at the projection plane, w is +inf
                    //
                    // 2 is a good recursion depth to appreciate the reflections without
                    // slowing things down too much
                    //
// window.ct= window.ct ? window.ct + 1 : 1; console.log(window.ct, c, x, y, w, h);
                    heap8[out_idx++] = trace_ray(0, 1, 0, x / w, y / w, 1, 1, w, 2);
                }
                heap8[out_idx++] = 255; // Opaque alpha
            }
        }

        heap8[inx_output]= out_idx - inx_output - 1;
        heap[inx_output]= out_idx - inx_output - 1;
        
        // console.log(heap8[inx_output], heap[inx_output]);
    };

    return {
        run: run
    };
}
