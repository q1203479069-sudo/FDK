function volume = fdk_reconstruction(proj, geo)
%FDK_RECONSTRUCTION Feldkamp-Davis-Kress (FDK) cone-beam CT reconstruction.
%   volume = fdk_reconstruction(proj, geo)
%
% Inputs
%   proj : [nu, nv, nBeta] projection data.
%          proj(:,:,k) is the detector image at angle beta(k).
%   geo  : struct with fields
%       .DSO        Distance source -> rotation center (mm)
%       .DSD        Distance source -> detector (mm)
%       .du         Detector pixel size along u (mm)
%       .dv         Detector pixel size along v (mm)
%       .nu         Number of detector pixels along u
%       .nv         Number of detector pixels along v
%       .nx         Reconstructed volume size along x
%       .ny         Reconstructed volume size along y
%       .nz         Reconstructed volume size along z
%       .dx         Voxel size along x (mm)
%       .dy         Voxel size along y (mm)
%       .dz         Voxel size along z (mm)
%       .beta       Projection angles (rad), length = nBeta
%       .filterType (optional) 'ram-lak' (default), 'shepp-logan', 'cosine'
%
% Output
%   volume : [nx, ny, nz] reconstructed 3D volume.
%
% Notes
%   1) Geometry assumes a circular scan with flat panel detector.
%   2) Coordinate system:
%      - Source rotates around z-axis.
%      - Detector u axis is tangent to orbit, v axis is z direction.

arguments
    proj (:,:,:) double
    geo struct
end

validate_geometry(geo, size(proj));

if ~isfield(geo, 'filterType')
    geo.filterType = 'ram-lak';
end

nu = geo.nu; nv = geo.nv; nBeta = numel(geo.beta);

% Detector coordinates centered at zero.
u = ((0:nu-1) - (nu-1)/2) * geo.du;
v = ((0:nv-1) - (nv-1)/2) * geo.dv;
[U, V] = ndgrid(u, v);

% ---- Step 1: Distance weighting ----
% Weight = DSD / sqrt(DSD^2 + u^2 + v^2)
W = geo.DSD ./ sqrt(geo.DSD^2 + U.^2 + V.^2);
weightedProj = proj .* W;

% ---- Step 2: 1D filtering along detector-u direction ----
filteredProj = zeros(size(weightedProj), 'like', weightedProj);
for k = 1:nBeta
    filteredProj(:,:,k) = filter_projection_u(weightedProj(:,:,k), geo.du, geo.filterType);
end

% ---- Step 3: Backprojection ----
volume = zeros(geo.nx, geo.ny, geo.nz, 'double');

x = ((0:geo.nx-1) - (geo.nx-1)/2) * geo.dx;
y = ((0:geo.ny-1) - (geo.ny-1)/2) * geo.dy;
z = ((0:geo.nz-1) - (geo.nz-1)/2) * geo.dz;
[Y, X] = meshgrid(y, x); % X/Y in image coordinates

for iz = 1:geo.nz
    z0 = z(iz);
    sliceAcc = zeros(geo.nx, geo.ny);

    for k = 1:nBeta
        beta = geo.beta(k);


        % Ray denominator: distance from voxel to source projected to detector normal.
        t = geo.DSO - X * cos(beta) - Y * sin(beta);

        % Avoid division by zero or negative values (outside valid cone).
        valid = t > 1e-6;
        if ~any(valid, 'all')
            continue;
        end

        % Detector coordinates (u,v) where voxel projects.
        uStar = geo.DSD * (-X * sin(beta) + Y * cos(beta)) ./ t;
        vStar = geo.DSD * (z0) ./ t;

        % Map physical coordinates to detector index space for interpolation.
        iu = (uStar / geo.du) + (nu+1)/2;
        iv = (vStar / geo.dv) + (nv+1)/2;

        % Bilinear interpolation from filtered projection.
        projk = filteredProj(:,:,k);
        val = interp2(1:nv, 1:nu, projk, iv, iu, 'linear', 0);

        % FDK backprojection factor (magnification correction).
        val = val .* (geo.DSO ./ t).^2;
        val(~valid) = 0;

        sliceAcc = sliceAcc + val;
    end

    % Discrete integration over beta.
    dBeta = mean(diff(unwrap(geo.beta)));
    volume(:,:,iz) = sliceAcc * dBeta;
end

end

function out = filter_projection_u(p, du, filterType)
%FILTER_PROJECTION_U Apply ramp-like filter along first dimension (u).

[nu, nv] = size(p);

nfft = 2^nextpow2(2*nu);
freq = (0:nfft-1)' / nfft; % normalized [0,1)
freq(freq > 0.5) = freq(freq > 0.5) - 1; % shift to [-0.5,0.5)
omega = abs(freq);

% Base ramp filter in frequency domain.
H = 2 * omega;

switch lower(filterType)
    case 'ram-lak'
        % no window
    case 'shepp-logan'
        f = omega;
        nz = f > 0;
        H(nz) = H(nz) .* sin(pi*f(nz)/(2*max(f(nz)))) ./ (pi*f(nz)/(2*max(f(nz))));
    case 'cosine'
        H = H .* cos(pi*omega/(2*max(omega)));
    otherwise
        error('Unknown filterType: %s', filterType);
end

P = fft(p, nfft, 1);
P = P .* H;
out = real(ifft(P, [], 1));
out = out(1:nu, :) / du;

end

function validate_geometry(geo, projSize)
required = {'DSO','DSD','du','dv','nu','nv','nx','ny','nz','dx','dy','dz','beta'};
for i = 1:numel(required)
    if ~isfield(geo, required{i})
        error('geo.%s is required.', required{i});
    end
end

if projSize(1) ~= geo.nu || projSize(2) ~= geo.nv
    error('Projection size mismatch: proj is [%d,%d,*], but geo expects [%d,%d,*].', ...
        projSize(1), projSize(2), geo.nu, geo.nv);
end

if numel(geo.beta) ~= projSize(3)
    error('Projection angle count mismatch: numel(beta)=%d, nProj=%d.', ...
        numel(geo.beta), projSize(3));
end

if geo.DSD <= geo.DSO
    warning('Usually DSD > DSO for flat-panel CBCT geometry.');
end

end
