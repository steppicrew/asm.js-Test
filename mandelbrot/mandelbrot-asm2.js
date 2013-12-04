

function MandelInner() {
    'use asm';

    function data( MaxIterations, c_re, c_im  ) {
        MaxIterations= MaxIterations|0;
        c_re = +c_re;
        c_im = +c_im;

        var n = 0;
        var Z_re= 0.0;
        var Z_im= 0.0;
        var Z_re2= 0.0;
        var Z_im2= 0.0;

        Z_re= +c_re;
        Z_im= +c_im;

        for ( n = 0; (n|0) < (MaxIterations|0); n= (n + 1)|0 ) {
            Z_re2 = Z_re * Z_re;
            Z_im2 = Z_im * Z_im;
            if ( Z_re2 + Z_im2 > +4 ) return n|0;

            Z_im = +2 * Z_re * Z_im + c_im;
            Z_re = Z_re2 - Z_im2 + c_re;
        }
        return 0;
    }

    return {
        data: data,
    }
}

var mandelInner= MandelInner();

L.MandelbrotSet = L.TileLayer.Canvas.extend({
    tileSize: 256,

    initialize: function (options) {
        L.Util.setOptions(this, options);
        this.drawTile = function (canvas, tilePoint) {
            var ctx = {
                canvas: canvas,
                tile: tilePoint
            };
            this._draw(ctx);
        };

        var profileTimer;
        this.on('loading', function() {
            profileTimer= new Date();
        })
        this.on('load', function() {
            document.getElementById('profiler').innerHTML= (new Date() - profileTimer) + ' ms<br />';
        })
    },

    _data: function(ctx) {
        var tileCount = 1 << this._map._zoom;
        var ReStart = -2.0;
        var ReDiff = 3.0;
        var MinRe = ReStart + ReDiff * ctx.tile.x / tileCount;
        var MaxRe = MinRe + ReDiff / tileCount;
        var ImStart = -1.2;
        var ImDiff = 2.4;
        var MinIm = ImStart + ImDiff * ctx.tile.y / tileCount;
        var MaxIm = MinIm + ImDiff / tileCount;
        var Re_factor = (MaxRe - MinRe) / (this.tileSize - 1);
        var Im_factor = (MaxIm - MinIm) / (this.tileSize - 1);
        var MaxIterations = 32;
        
        var data = [];
        for (var y = 0, i = 0; y < this.tileSize; ++y) {
            var c_im = MinIm + y * Im_factor;
            for (var x = 0; x < this.tileSize; ++x) {
                var c_re = MinRe + x * Re_factor;
                var n= mandelInner.data(MaxIterations, c_re, c_im);

                if (n < MaxIterations / 2) {
                    data[i++] = 255 / (MaxIterations / 2) * n;
                    data[i++] = data[i++] = 0;
                } else {
                    data[i++] = 255;
                    data[i++] = data[i++] = (n - MaxIterations / 2) * 255 / (MaxIterations / 2);
                }
                data[i++] = 255;
            }
        }        
        return data;
    },
    
    _draw: function (ctx) {
        var data = this._data(ctx);
        var g = ctx.canvas.getContext('2d');
        var d = g.createImageData(this.tileSize, this.tileSize);
        for (var i = 0, n = this.tileSize * this.tileSize * 4; i < n; i++)  {
            d.data[i] = data[i];
        }
        g.putImageData(d, 0, 0);
    }    
});

var attr = 'Adapted from a <a target="_blank" href="http://polymaps.org/ex/mandelbrot.html">polymaps sample</a>'
var map = L.map('map');
map.attributionControl.addAttribution(attr);
map.addLayer(new L.MandelbrotSet());
map.setView(new L.LatLng(0, 0), 2);
