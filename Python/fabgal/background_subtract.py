import numpy as np
from scipy.ndimage import uniform_filter, grey_opening
from skimage.transform import resize
from skimage.color import rgb2hsv, hsv2rgb

class BackgroundSubtracter:
    """
    Reference Python implementation of ImageJ's BackgroundSubtracter (v5.0).
    
    Geometric Assumption:
    The "Rolling Ball" algorithm assumes that the Z-axis (intensity) and X/Y-axes (spatial) 
    share a common unit scale (1.0 intensity ≈ 1.0 pixel). 
    - For integer images (counts), this is usually physically valid.
    - For normalized float images [0.0, 1.0], a radius of 50px creates a ball 50.0 units deep,
      which is astronomically larger than the data range. In these cases, use 'ball_z_factor'.
    """

    def run(self, image: np.ndarray, radius: float, 
            light_background: bool = True, 
            use_paraboloid: bool = False, 
            do_presmooth: bool = True,
            rgb_mode: str = 'independent',
            output_type: str = 'difference',
            float_pivot: float = None,
            ball_z_factor: float = 1.0,
            nan_policy: str = 'propagate'):
        """
        Args:
            image (np.ndarray): Input image.
            radius (float): Rolling ball radius in pixels.
            light_background (bool): True for brightfield, False for fluorescence.
            use_paraboloid (bool): Use 2D parabolic kernel (smoother).
            do_presmooth (bool): Apply 3x3 mean filter before background estimation.
            rgb_mode (str): 'independent' or 'hsb'.
            output_type (str): 'difference' (signal), 'corrected' (flattened), 'background'.
            float_pivot (float): Inversion pivot for float images. Default: data max.
            ball_z_factor (float): Scale factor for intensity vs spatial units. 
                                   Example: If 1.0 intensity unit = 100 pixels, set to 0.01.
                                   Useful for normalized float images [0,1].
            nan_policy (str): 'propagate' (keep NaNs), 'omit' (fill with local min/zero).
        """
        image = np.asanyarray(image)
        
        # 0. Guard Clauses
        if radius < 1:
            return image.copy()
        
        # 1. RGB Handling
        if image.ndim == 3 and image.shape[-1] in [3, 4]:
            
            has_alpha = (image.shape[-1] == 4)
            rgb = image[..., :3]
            alpha = image[..., 3] if has_alpha else None

            if rgb_mode == 'hsb':
                # HSB Mode: Recurse on V channel
                rgb_norm, scale_factor = self._normalize_for_color_transform(rgb)
                rgb_norm = np.clip(rgb_norm, 0, 1)
                hsv = rgb2hsv(rgb_norm)
                
                # Pass explicit z-factor logic to recursion
                hsv[..., 2] = self.run(hsv[..., 2], radius, light_background, 
                                       use_paraboloid, do_presmooth, 
                                       output_type=output_type,
                                       float_pivot=1.0,
                                       ball_z_factor=ball_z_factor,
                                       nan_policy=nan_policy)
                
                rgb_out_norm = hsv2rgb(hsv)
                rgb_out = self._denormalize_color_transform(rgb_out_norm, scale_factor, image.dtype)
                if has_alpha: return np.dstack((rgb_out, alpha))
                return rgb_out
            else:
                # Independent Mode
                processed = [self.run(rgb[..., i], radius, light_background, 
                                      use_paraboloid, do_presmooth, rgb_mode, 
                                      output_type, float_pivot, ball_z_factor, nan_policy) 
                             for i in range(3)]
                rgb_out = np.stack(processed, axis=-1)
                if has_alpha: return np.dstack((rgb_out, alpha))
                return rgb_out

        # 2. Data Preparation
        original_dtype = image.dtype
        img_float = image.astype(np.float64, copy=False)

        # 2a. NaN Policy
        if nan_policy == 'omit':
            # Replace NaNs with a "neutral" value (min for background logic)
            # This prevents NaNs from poisoning the morphological window
            mask_nan = np.isnan(img_float)
            if np.any(mask_nan):
                fill_val = np.nanmin(img_float) if not light_background else np.nanmax(img_float)
                img_float[mask_nan] = fill_val

        # 3. Domain Inversion & Z-Scaling
        # Determine pivot
        if np.issubdtype(original_dtype, np.integer):
            pivot = np.iinfo(original_dtype).max
        else:
            pivot = float_pivot if float_pivot is not None else (np.nanmax(img_float) or 1.0)

        # Invert topology
        if light_background:
            img_float = pivot - img_float
        
        # Apply Z-Factor (Geometric Scaling)
        # We scale the "height" of the image features relative to the ball radius.
        if ball_z_factor != 1.0:
            img_float *= ball_z_factor
            pivot *= ball_z_factor

        # 4. Pre-smoothing (Reflect mode)
        if do_presmooth:
            img_float = uniform_filter(img_float, size=3, mode='reflect')

        # 5. Adaptive Downscaling
        shrink_factor = 1
        if radius > 100: shrink_factor = 8
        elif radius > 30: shrink_factor = 4
        elif radius > 10: shrink_factor = 2
        
        # 6. Edge Handling (Pad)
        orig_shape = img_float.shape
        pad_y = (shrink_factor - (orig_shape[0] % shrink_factor)) % shrink_factor
        pad_x = (shrink_factor - (orig_shape[1] % shrink_factor)) % shrink_factor
        if pad_y > 0 or pad_x > 0:
            img_float = np.pad(img_float, ((0, pad_y), (0, pad_x)), mode='reflect')

        # 7. Create Kernel
        # Ensure kernel radius is at least 1.0 to maintain structural integrity
        scaled_radius = max(radius / shrink_factor, 1.0)
        kernel = self._create_kernel(scaled_radius, use_paraboloid)

        # 8. Downscale (Min/Erosion)
        if shrink_factor > 1:
            h, w = img_float.shape
            sh, sw = h // shrink_factor, w // shrink_factor
            reshaped = img_float.reshape(sh, shrink_factor, sw, shrink_factor)
            small_img = reshaped.min(axis=(1, 3))
        else:
            small_img = img_float

        # 9. Morphology (Background Hull)
        background_small = grey_opening(small_img, structure=kernel, mode='reflect')

        # 10. Upscale
        if shrink_factor > 1:
            background = resize(background_small, img_float.shape, order=1, 
                                preserve_range=True, anti_aliasing=False)
        else:
            background = background_small

        # 11. Crop & Subtract
        if pad_y > 0 or pad_x > 0:
            background = background[:orig_shape[0], :orig_shape[1]]
            img_float = img_float[:orig_shape[0], :orig_shape[1]]
        
        signal = img_float - background
        signal = np.clip(signal, 0, None)
        
        # 12. Output Logic (and De-scaling Z)
        if ball_z_factor != 1.0:
            signal /= ball_z_factor
            background /= ball_z_factor
            pivot /= ball_z_factor
            # img_float (smoothed input) was modified in place, conceptually we revert scale
            
        if output_type == 'difference':
            result = signal
        elif output_type == 'background':
            result = (pivot - background) if light_background else background
        elif output_type == 'corrected':
            result = (pivot - signal) if light_background else signal
        else:
            raise ValueError(f"Unknown output_type: {output_type}")

        return self._restore_dtype(result, original_dtype)

    def _create_kernel(self, radius, use_paraboloid):
        radius_int = int(np.ceil(radius))
        x, y = np.meshgrid(np.arange(-radius_int, radius_int + 1),
                           np.arange(-radius_int, radius_int + 1))
        dist_sq = x**2 + y**2

        if use_paraboloid:
            z = dist_sq / (2 * radius)
            kernel = -z
        else:
            mask = dist_sq <= radius**2
            z = np.zeros_like(dist_sq, dtype=np.float64)
            z[mask] = np.sqrt(radius**2 - dist_sq[mask])
            kernel = z.copy()
            kernel[~mask] = -np.inf
            
        kernel -= kernel.max()
        return kernel

    def _normalize_for_color_transform(self, rgb):
        orig_dtype = rgb.dtype
        rgb_f = rgb.astype(np.float64, copy=False)
        if np.issubdtype(orig_dtype, np.integer):
            scale = np.iinfo(orig_dtype).max
        else:
            scale = np.nanmax(rgb_f)
            if scale == 0: scale = 1.0
        return rgb_f / scale, scale

    def _denormalize_color_transform(self, rgb_norm, scale, dtype):
        rgb = rgb_norm * scale
        return self._restore_dtype(rgb, dtype)

    def _restore_dtype(self, image, dtype):
        if np.issubdtype(dtype, np.integer):
            info = np.iinfo(dtype)
            image = np.clip(image, info.min, info.max)
            return np.rint(image).astype(dtype)
        else:
            return image.astype(dtype, copy=False)