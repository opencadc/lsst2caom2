# # Butler-2-CAOM2
import caom2
import datetime
from lsst.daf.butler import Butler, DatasetRef
import os
import requests
from upath import UPath
import urllib
from astropy.coordinates import EarthLocation

# A bunch of static variables that define our problem.  These are set to work
# for the LSSTComCam DP1 release.
READ_GROUPS = caom2.caom_util.URISet(
    "ivo://cadc.nrc.ca/gms?LSST_CADC", "ivo://cadc.nrc.ca/gms?LSST_Can"
)
location = EarthLocation.of_site("Rubin")
TELESCOPE = caom2.Telescope(
    name="8.4-meter Simonyi Survey Telescope",
    geo_location_x=location.geocentric[0].to("m").value,
    geo_location_y=location.geocentric[1].to("m").value,
    geo_location_z=location.geocentric[2].to("m").value,
)
INSTRUMENT = "LSSTComCam"
OBSTYPE = "OBJECT"
CAOM_COLLECTION = "LSST.DP1"
LSST_PROPOSAL_TITLE = "Legacy Survey of Space and Time"
LSST_PROJECT = "LSST.DP1"
LSST_PROPOSAL_ID = "DP1"
BUTLER = "dp1"
LSST_COLLECTIONS = "LSSTComCam/DP1"
PROPOSAL = "LSST"
SKYMAP = "lsst_cells_v1"
CONTENT_TYPE = "application/fits"
OBSERVABLE = caom2.Observable("phot.flux.density")
OBSERVATION_ID_KEYS = ["skymap", "tract", "patch"]
DATASET_TYPE = "deep_coadd"
ALGORITHM = f"lsst.{DATASET_TYPE}.DRP.DP1.DM-51335"
SI_PATH_ROOT = "/raven/files/"
META_RELEASE = datetime.datetime.now(datetime.UTC)
PROVENANCE_NAME = "LSST Science Pipeline"
PRODUCER = "Rubin Observatory"
PROVENANCE_REFERENCE = "https://dp1.lsst.io/index.html"
# From the SVO which took it
# from: https://github.com/lsst/throughputs/tree/main/baseline
# values are in metres
BANDPASS_TOTAL_THROUGHPUT = {
    "u": [3206.34E-10, 4081.51E-10],
    "g": [3876.02E-10, 5665.33E-10],
    "r": [5377.19E-10, 7055.16E-10],
    "i": [6765.77E-10, 8325.05E-10],
    "z": [8035.39E-10, 9375.47E-10],
    "y": [9089.07E-10, 10915.01E-10],
}
PIXEL_SCALE = 0.2 / 3600.0  # pixel scale in degrees.
COADD_DIMENSION = caom2.Dimension2D(3400, 3400)


def get_content_length(butler: Butler, dataset_ref: DatasetRef) -> int:
    """
    To label each artifact with its content-length the easiest approach appears
    be to use the SI header.
    """
    response = requests.head(
        butler.getURI(dataset_ref).geturl(),
        headers={"Authorization": f"Bearer {os.getenv('CADC_TOKEN')}"},
    )
    return int(response.headers["content-length"])


def get_parts(butler: Butler,
              dataset_ref: DatasetRef) -> caom2.TypedOrderedDict:
    """
    return an ordered list of parts for a standard LSST coadd images
    The coadd components of the LSST DRP are a fixed structure. Here we define
    3 parts for the FITS extension which are the standard image/mask/variance
    group of the LSST DRP coadd.  The other extensions contain metadata and
    are not enumerated as accessible parts.
    """
    if dataset_ref.datasetType.name != DATASET_TYPE:
        return None
    _image = butler.get(dataset_ref)
    wcs = _image.getWcs()
    naxis1, naxis2 = _image.getDimensions()
    [[cd11, cd12], [cd21, cd22]] = wcs.getCdMatrix()
    crpix1, crpix2 = wcs.getPixelOrigin()
    crval1, crval2 = wcs.getSkyOrigin()
    # pixel_scale = wcs.getPixelScale()
    cdmatrix = caom2.CoordFunction2D(
        dimension=caom2.Dimension2D(naxis1, naxis2),
        ref_coord=caom2.Coord2D(
            caom2.RefCoord(crpix1, crval1.asDegrees()),
            caom2.RefCoord(crpix2, crval2.asDegrees()),
        ),
        cd11=cd11,
        cd12=cd12,
        cd21=cd21,
        cd22=cd22,
    )
    position = caom2.SpatialWCS(
        axis=caom2.wcs.CoordAxis2D(
            axis1=caom2.Axis(ctype="RA-TAN", cunit="deg"),
            axis2=caom2.Axis(ctype="DEC-TAN", cunit="deg"),
            function=cdmatrix,
        ),
        coordsys="ICRS",
        equinox=2000.0,
    )
    product_types = {
        "IMAGE": caom2.ProductType.THIS,
        "MASK": caom2.ProductType.AUXILIARY,
        "VARIANCE": caom2.ProductType.AUXILIARY,
    }
    band = dataset_ref.dataId['band']
    lwave = BANDPASS_TOTAL_THROUGHPUT[band][0]
    uwave = BANDPASS_TOTAL_THROUGHPUT[band][1]
    delta = uwave-lwave
    eff = (uwave+lwave)/2.0
    energy = caom2.SpectralWCS(axis=caom2.wcs.CoordAxis1D(
        axis=caom2.Axis(ctype='WAVE', cunit='m'),
        error=caom2.CoordError(1E-9, 1E-9),
        function=caom2.CoordFunction1D(naxis=1,
                                       delta=delta,
                                       ref_coord=caom2.RefCoord(1.0, eff))),
        specsys='TOPOCENT',
        ssysobs='TOPOCENT',
        ssyssrc='TOPOCENT',
        resolving_power=eff/delta,
        bandpass_name=band)
    parts = caom2.TypedOrderedDict(caom2.Part)
    for part_name in product_types:
        product_type = product_types[part_name]
        chunk = caom2.Chunk(
            product_type=product_type,
            naxis=2,
            position_axis_1=1,
            position_axis_2=2,
            position=position,
            energy_axis=None,
            energy=energy
        )
        part = caom2.Part(
            part_name,
            product_type=product_type,
            chunks=caom2.TypedList(caom2.Chunk, chunk),
        )
        parts.add(part)
    return parts


def get_artifacts(butler: Butler,
                  dataset_type: str,
                  instrument: str,
                  data_id: {}) -> [caom2.TypedOrderedDict]:
    """
    Given a patch of the teselation of the sky retrieve a list of Artifact
    instances.
    butler: the LSST Data Butler where this data is registered.
    dataset_type: the name of coadd (e.g. deep_coadd).
    instrument:  the instrument that aquired this data.
    data_id:  dimensions that identify the desired coadd for this instrument.
    get_artifacts expects the butler will have 'coadd', 'coadd_n_image' and
    'coadd_background' products.
    """
    product_types = {
        f"{dataset_type}": caom2.ProductType.THIS,
        f"{dataset_type}_n_image": caom2.ProductType.AUXILIARY,
        f"{dataset_type}_background": caom2.ProductType.AUXILIARY,
    }
    artifacts = caom2.TypedOrderedDict(caom2.Artifact)
    content_type = CONTENT_TYPE
    dataset_refs = butler.query_all_datasets(data_id=data_id)
    for dataset_ref in dataset_refs:
        if dataset_ref.datasetType.name not in product_types:
            continue
        uri = urllib.parse.urlsplit(butler.getURI(dataset_ref).geturl()
                                    ).path.replace(SI_PATH_ROOT, "")
        artifacts.add(
            caom2.Artifact(
                uri=uri,
                product_type=product_types[dataset_ref.datasetType.name],
                release_type=caom2.ReleaseType.DATA,
                content_type=content_type,
                content_length=get_content_length(butler, dataset_ref),
                parts=get_parts(butler, dataset_ref),
            )
        )
    return artifacts


def get_plane_position(butler: Butler, data_id: dict) -> caom2.Position:
    """
    Get the position region of this data Id.
    """
    patches = butler.query_dimension_records("patch", data_id=data_id)
    assert len(patches) == 1
    ivoa_pos = patches[0].region.to_ivoa_pos()
    shapes = {
        "POLYGON": caom2.shape.Polygon,
        "CIRCLE": caom2.shape.Circle,
        "BOX": caom2.shape.Box,
    }
    shape = shapes[ivoa_pos.split()[0]](ivoa_pos.split()[1:])
    return caom2.Position(
        bounds=shape, dimension=COADD_DIMENSION, sample_size=PIXEL_SCALE
    )


def get_plane_energy(butler: Butler, data_id: dict) -> caom2.Energy:
    """
    Get the energy bounds given the data_id
    """
    interval = BANDPASS_TOTAL_THROUGHPUT[data_id["band"]]
    R = ((interval[0] + interval[1]) / 2) / (interval[1] - interval[0])
    interval = caom2.Interval(lower=interval[0],
                              upper=interval[1])
    return caom2.Energy(
        bandpass_name=data_id["band"],
        bounds=interval,
        resolving_power=R,
        em_band=caom2.EnergyBand.OPTICAL,
    )


def get_planes(
    butler: Butler, dataset_type: str, instrument: str, data_id: dict
) -> [caom2.Plane]:
    """
    Build plane objects for each dataset_ref in the list and return a set of
    caom2.Plane object.

    get_planes build one plane for each 'dataset_type' and to put all the the
    artifacts with the same data_id int that plane.
    
    The list of dataset_refs all have the same data_id and its expected that
    band is NOT in data_id and all have datasetType.name==DATASET_TYPE

    A list of planes is returned.
    """
    planes = caom2.TypedOrderedDict(caom2.Plane)
    for dataset_ref in butler.query_datasets(dataset_type, 
                                             instrument=instrument, 
                                             data_id=data_id):
        # The dataset_ref will have the same data_id as above augemented with
        # the 'band' value. All datasets with this data_id (restricted by band)
        # are in the same artifact
        band = dataset_ref.dataId['band']
        product_id = f"{dataset_type}-{band}"
        artifacts = get_artifacts(butler,
                                  dataset_type,
                                  instrument,
                                  dataset_ref.dataId)
        data_read_groups = READ_GROUPS
        data_product_type = caom2.DataProductType.IMAGE
        calibration_level = caom2.CalibrationLevel.PRODUCT
        run = UPath(dataset_ref.run).parts
        provenance = caom2.Provenance(
            name=PROVENANCE_NAME,
            version=run[-1],
            project=run[2],
            producer=PRODUCER,
            run_id=dataset_ref.run,
            reference=PROVENANCE_REFERENCE,
        )
        plane = caom2.Plane(
            product_id=product_id,
            artifacts=artifacts,
            data_read_groups=data_read_groups,
            data_product_type=data_product_type,
            calibration_level=calibration_level,
            provenance=provenance,
            observable=OBSERVABLE,
        )
        # plane.position = get_plane_position(butler, dataset_ref.dataId)
        # plane.energy = get_plane_energy(butler, dataset_ref.dataId)
        planes[product_id] = plane
    return planes


def get_observation(butler: Butler,
                    dataset_type: str,
                    instrument: str,
                    data_id: dict) -> caom2.Observation:
    """
    butler: LSST Data Butler object
    dataset_type: dataset that will be the primary artifact for each plane
    instrument: the name of the instrument defines the schema
    data_id: a dictionary defining the sky patch that is this observation.

    Build the caom2 record for this skymap/tract/patch set.
    skymap/tract/patch sets the 'observation_id' of a stack.
    """
    observation_id = "-".join([f"{data_id[x]}" for x in data_id])
    planes = get_planes(butler, dataset_type, instrument, data_id)
    proposal = caom2.Proposal(
        id=PROPOSAL, project=LSST_PROJECT, title=LSST_PROPOSAL_TITLE
    )
    instrument = caom2.Instrument(INSTRUMENT)
    return caom2.Observation(
        collection=CAOM_COLLECTION,
        observation_id=observation_id,
        algorithm=ALGORITHM,
        intent=caom2.ObservationIntentType.SCIENCE,
        type=OBSTYPE,
        proposal=proposal,
        telescope=TELESCOPE,
        instrument=instrument,
        meta_release=META_RELEASE,
        planes=planes,
    )


def get_observation_ids(butler: Butler,
                        dataset_type: str,
                        skymap: str,
                        keys: list):
    dataset_refs = butler.query_datasets(dataset_type, skymap=skymap)
    data_ids = []
    for dataset_ref in dataset_refs:
        data_id = {x: dataset_ref.dataId[x] for x in keys}
        if data_id not in data_ids:
            data_ids.append(data_id)
    return data_ids


if __name__ == "__main__":
    butler = Butler(BUTLER, collections=LSST_COLLECTIONS)
    for data_id in get_observation_ids(butler,
                                       DATASET_TYPE,
                                       SKYMAP,
                                       OBSERVATION_ID_KEYS):
        record = get_observation(butler, DATASET_TYPE, INSTRUMENT, data_id)
        break
    print(record)
